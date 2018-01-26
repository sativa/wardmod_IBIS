c hourlymet.f
c Created by Jason Patton 2010-03-25
c
c Holding place for subprograms dealing with using hourly met data
c
c Quick notes:
c   Hourly files are opened/read/closed using UNIT = 13
c   Files are named by year
c   The order of variables in files must be:
c     year month day hour temp rh wind precip solar pres

c ---------------------------
      FUNCTION getUnit(varName)
c ---------------------------
c Created by Jason Patton 2010-03-24
c
c Lets user define units for variables for 
c hourly met data.  These units are then
c retrieved using getUnit function.
c
c Elevation is also defined
c
c     Output var
      CHARACTER*5 getUnit

c     Input var
      CHARACTER*4 varName

c     Processing var
      CHARACTER*5 units(8)

c Temperature (C, K, or F)
      units(1) = 'F'

c Humidity
c     dewpC = Dewpoint C
c     dewpK = Dewpoint K
c     dewpF = Dewpoint F
c     rhpct = RH Percent (0% - 100%)
c     rhfrc = RH Fraction (0.0 - 1.0)
c     mixgg = Mixing Ratio (g/g)
c     mxgkg = Mixing Ratio (g/kg)
      units(2) = 'rhpct'

c Wind (ms = m/s, kmph = km/h, mph = mi/h, or kts = nm/h)
      units(3) = 'mph'

c Precip (mm, cm, or in)
      units(4) = 'mm'

c Solar Radiation (Wm2 = W m-2, kcal = kcal m-2 h-1)
      units(5) = 'kcal'

c Pressure (mb, kPa, in)
      units(6) = 'in'

c Carbon Dioxide (ppm, mgm3 = mg/m3)
c      units(7) = 'ppm'

c Spit back units as function
      IF(varName == 'Temp') getUnit = units(1)
      IF(varName == 'Humi') getUnit = units(2)
      IF(varName == 'Wind') getUnit = units(3)
      IF(varName == 'Prec') getUnit = units(4)
      IF(varName == 'Sola') getUnit = units(5)
      IF(varName == 'Pres') getUnit = units(6)

      RETURN

c ----------------------------------------------------
      END FUNCTION getUnit
c ----------------------------------------------------
c
c ----------------------------------------------------
      SUBROUTINE openHourlyMet (thisYear, seed, doSpinup)
c ----------------------------------------------------
c Created by Jason Patton 2010-03-25
c
c This should run at the beginning of every year
c
c 1. Opens hourly data files named by year with UNIT = 13
c 2. Checks for open errors and stops program if error found

c     Force explicitly defined variables
      IMPLICIT NONE

c     Constants
      INTEGER availYears(21) ! Change these based on years available
      DATA availYears / 1989, 1990, 1991, 1992, 1993, 1994, 1995,
     >                  1996, 1997, 1998, 1999, 2000, 2001, 2002,
     >                  2003, 2004, 2005, 2006, 2007, 2008, 2009 /

c     Incoming vars
      INTEGER thisYear
      INTEGER seed
      INTEGER doSpinup

c     Processing vars
      LOGICAL isLeap, is_leap
      INTEGER i, oer, metYear, randYear
      CHARACTER*4 charfile
      CHARACTER*4 checkYear

c Set curent year
      metYear = thisYear
      randYear = 0

c Check for leapyear
      isLeap = IS_LEAP(metYear)

c If this is a spinup case, get a random year of data
      IF(doSpinup .EQ. 1) THEN
         DO
            CALL RANDOM_CHOICE(randYear, availYears, 21, seed)
c Loop until we get a match for leap year/not leap year
            IF(                                                ! If...
     >        ((      isLeap) .AND. (      IS_LEAP(randYear))) ! both years are leap years...
     >                        .OR.                             !         ...or...
     >        ((.NOT. isLeap) .AND. (.NOT. IS_LEAP(randYear))) ! both years are not leap years...
     >      ) EXIT                                             ! the program may continue.
         ENDDO
         WRITE(*,*) 'Spinup run, using data from year (', randYear,')'
      ENDIF

c Tell user that they're using hourly data this year
      WRITE(*,*) 'Using hourly data for year (',metYear,')'

c Put year into string for filename
      IF(randYear .NE. 0) THEN
         WRITE(charfile,'(I4)') randYear
      ELSE
         WRITE(charfile,'(I4)') metYear
      ENDIF

c Open file
      OPEN (UNIT = 13, FILE = 'input/hourly/'//charfile, IOSTAT=oer)

c Check for error
      IF (oer .NE. 0) THEN
         WRITE(*,*) 'Error: Hourly data file unopenable for:'
         WRITE(*,*) '  ', metYear
         STOP
      ENDIF

c Check for header, if not there, rewind file
      READ(13,*) checkYear
      IF(checkYear .EQ. charfile) REWIND(13)

      RETURN

c ----------------------------------------------------
      END SUBROUTINE openHourlyMet
c ----------------------------------------------------
c
c ----------------------------------------------------
      SUBROUTINE readHourlyMet (thisYear, thisMonth, thisDay,
     .     thisStep, thisDOY, doirrigate, irrlength, irrstart,
     .     irrend)
c ----------------------------------------------------
c
c Created by Jason Patton 2010-03-25
c
c Call this at the beginning of each time step
c
c 1. Reads next line in hourly data file for this year (UNIT = 13).
c 2. Checks read data for date/time consistency
c 3. Converts data to IBIS-friendly units
c 4. Pushes data to IBIS's LSM vars

c     Get needed modules
      USE comgrid
      USE compar
      USE comwork
      USE comatm
      USE comhyd
      USE comsum
      USE combcs
      USE comcrop
      USE comnitr

c     Force explicitly defined variables
      IMPLICIT NONE

c     Incoming vars
      INTEGER thisYear, thisMonth, thisDay, thisStep,
     >     thisDOY, doirrigate
      real:: irrlength,irrstart, irrend

c     Read vars
      INTEGER metYear, metMonth, metDay, metHour
      REAL metTemp, metHum, metWind, metPrec,
     >     metSolar, metPres 
c     >     , metElev ! Uncomment if reading Elevation
c     >     , metCO2  ! Uncomment if reading CO2

c     Function vars
      CHARACTER*5 getUnit

c     Processing vars
      INTEGER ier,  ! Input status 
     >     i        ! General counter
      REAL metSH,   ! Computed specific humidity 
     >     metVP,   ! Computed vapor pressure
     >     minTemp, ! Computed minimum temp
     >     maxTemp, ! Computed maximum temp
     >     sumTemp, ! Sum of all temps for day
     >     avgTemp, ! Average temp for day
     >     sumSH,   ! Sum of all specific humidity for day
     >     avgSH,   ! Average specific humidity for day
     >     sumPrec, ! Sum of all precip for day
     >     sumWind, ! Sum of all wind speeds for day
     >     avgWind  ! Average wind speed for day
c     >     , sumYearPrecip, sumYearRain, sumYearSnow ! Testing


c     Save summation variables
      SAVE sumTemp, sumSH, sumPrec, sumWind
c     >     , sumYearPrecip ! Testing

c     Include functions for humidity calculations
      INCLUDE 'comsat.h'

c
c <BEGIN Read from file>
c
c Current file structure is:
c Year, Month, Day, Hour, Temp, Hum, Wind, Prec, Solar, Pres

      READ(13, *, IOSTAT = ier)
     >     metYear, metMonth, metDay, metHour,
     >     metTemp, metHum, metWind, metPrec,
     >     metSolar, metPres 
c     >     , metCO2 ! Uncomment if reading CO2

c Check for read error
      IF (ier .NE. 0) THEN
         WRITE(*,*) 'Error reading file for:'
         WRITE(*,*) thisYear
         WRITE(*,*) ' at:'
         WRITE(*,*) thisYear, thisMonth, thisDay, thisStep
         STOP
      ENDIF

c Check for date/time consistency
      IF ( (thisMonth .NE. metMonth) .OR.
     >     (thisDay   .NE. metDay)   .OR.
     >     (thisStep  .NE. (metHour + 1)) )
     >     THEN
         WRITE(*,*) 'Date mismatch at:'
         WRITE(*,*) '  IBIS Year', thisYear, 'Data Year', metYear
         WRITE(*,*) '  IBIS Month', thisMonth, 'Data Month', metMonth
         WRITE(*,*) '  IBIS Day', thisDay, 'Data Day', metDay
         WRITE(*,*) '  IBIS Step', thisStep, 'Data Hour', metHour
         STOP
      ENDIF

c
c <BEGIN Convert data to IBIS-friendly units>
c
c First, we'll grab elevation if needed (it's hidden away in a string)
c Should find a better way to do this later
c      READ(getUnit('Elev'), '(F)') metElev

c Temperature needs to be in K
      IF(getUnit('Temp') == 'F') THEN ! F to K
         metTemp = 273.16 + ((metTemp - 32.0)*(5.0/9.0))
      ENDIF

c Pressure needs to be in Pa
c Though this seems out of order, pressure is needed for specific humidity calc
      IF(getUnit('Pres') == 'in') THEN ! inHg to Pa
         metPres = 3386.39 * metPres
      ENDIF

c Humidity needs to be converted to specific humidity and vapor pres
      IF(getUnit('Humi') == 'rhpct') THEN ! % to *
c     Set RH to 99% if > 99%
         IF(metHum .GT. 99.0) metHum = 99.0
         metSH = QSAT((metHum/100.0) * ESAT(metTemp), metPres)
         metVP = (metHum/100.0) * ESAT(metTemp)
      ENDIF

c Precipitation needs to be in mm
c ***NOTE*** This assumes that precip data is FOR next hour, not FROM past hour
      IF(getUnit('Prec') == 'mm') THEN ! mm to mm
         metPrec = metPrec
      ENDIF

c Wind speed needs to be in m s-1 and between 2.5 and 25.0
      IF(getUnit('Wind') == 'mph') THEN ! mph to m s-1
         metWind = 0.44704 * metWind
      ENDIF
      metWind = MAX(2.5, MIN(25.0, metWind)) ! 2.5 < ua < 25.0

c Solar radiation needs to be in W m-2
      IF(getUnit('Sola') == 'kcal') THEN ! kcal hr-1 to W
         metSolar = 1.162 * metSolar
      ENDIF

c
c <BEGIN Push to IBIS weather vars>
c
      CALL CONST(ta, npoi, metTemp)
      CALL CONST(rh, npoi, metHum)
      CALL CONST(qa, npoi, metSH)
c      CALL CONST(ea, npoi, metVP)
      CALL CONST(ua, npoi, metWind)
      CALL CONST(precip, npoi, metPrec)
      CALL CONST(cloud, npoi, metSolar) ! weather.f subroutine diurnalmet uses cloud as insolation
      CALL CONST(rads, npoi, metSolar)
      CALL CONST(psurf, npoi, metPres)

c     This subroutine will calculate everything we don't have:
c       solar transmission, direct & diffuse
c       infrared transmission
c       snow/rain rates
c       irrigation (if needed)
      CALL hourlyMetCalc(thisStep, thisDOY, 
     >     doirrigate, irrlength, irrstart, irrend)

c
c <BEGIN Calculate daily vars (tmax, tmin, tavg)>
c
c First hour of day, reset min and max temps and sums
      IF(metHour .EQ. 0) THEN
         CALL CONST(tmin, npoi, metTemp)
         CALL CONST(tmax, npoi, metTemp)
         sumTemp = 0.0
         sumSH = 0.0
         sumPrec = 0.0
         sumWind = 0.0
      ENDIF

c Compare tmin/tmax to current value and set
      minTemp = MIN(tmin(1), metTemp)
      maxTemp = MAX(tmax(1), metTemp)
      CALL CONST(tmin, npoi, minTemp)
      CALL CONST(tmax, npoi, maxTemp)

c Sum up variables
      sumTemp = sumTemp + metTemp
      sumSH = sumSH + metSH
      sumPrec = sumPrec + metPrec
      sumWind = sumWind + metWind



c Last hour of day, compute averages and call dailymet
      IF(metHour .EQ. 23) THEN
         avgTemp = sumTemp/(24.0)             ! Average temp from hourly temps
c         avgTemp = (minTemp + maxTemp) / 2.0 ! Average temp from min/max
         avgSH = sumSH/(24.0)
         avgWind = sumWind/(24.0)
         CALL CONST(td, npoi, avgTemp)
         CALL CONST(qd, npoi, avgSH)
         CALL CONST(precip, npoi, sumPrec)
         CALL CONST(ud, npoi, avgWind)
         CALL dailyMetCalc(thisMonth, thisDay, 0)
      ENDIF

      RETURN

c ----------------------------------------------------
      END SUBROUTINE readHourlyMet
c ----------------------------------------------------
c
c ----------------------------------------------------
      SUBROUTINE closeHourlyMet ()
c ----------------------------------------------------
c Created by Jason Patton 2010-03-24
c
c This should run at end of each year
c
c 1. Closes currently open hourly data file (UNIT = 13).

c     Force explicitly defined variables
      IMPLICIT NONE

c Close file
      CLOSE(13)

      RETURN

c ----------------------------------------------------
      END SUBROUTINE closeHourlyMet
c ----------------------------------------------------

c Since the input files can't retrieve all of the necessary variables that
c the LSX portion of the model needs to run, additional weather functions 
c are called by the above subroutines.  These are derived (i.e. stolen)
c from the weather.f functions called when using climate data.

c ---------------------------------------------------------------------
      SUBROUTINE hourlyMetCalc (istep, jday, doirrigate, 
     >     ilens, starti, endi)
c ---------------------------------------------------------------------
c
c uses:
c
      USE comgrid
      USE compar
      USE comatm
      USE comveg
      USE comwork
      USE comhour
      USE comcrop
      USE comnitr
      use combcs,only: cmask
c
      IMPLICIT NONE
c
c Arguments
c
      INTEGER istep,            ! current time step (between 1 and niter)
     >     jday,                ! current day of year
     >     doirrigate           ! irrigation switch
c
      REAL
     >     ilens,               ! length of irrigation events
     >     starti,              ! start time of irrigiation events
     >     endi,                ! end time of irrigation events
     >     truecloud            ! cloudiness computed by atm transmission
c
      INTEGER i,                ! loop indice
     >     jj,                  ! latitude
     >     ib                   ! waveband number 1= visible, 2= near-IR
c
      REAL time,                ! time in seconds since midnight
     >     rtime,               ! time in hours
     >     orbit,               ! earth's orbital angle (around the sun) in radians
     >     xdecl,               ! solar declination angle
     >     sw,                  ! effective solar constant
     >     xlat,                ! latitude in radians
     >     fdiffuse,            ! fraction of indirect (diffuse) solar radiation
     >     wfrac,               ! fraction of energy in each waveband
     >     emb,
     >     ea,
     >     ec,
     >     dtair,
     >     dtcloud
c
      INCLUDE 'comsat.h'
c
c Externals
c
      REAL get_time             ! from utilities.f
      REAL calc_orbit, calc_xdecl, calc_sw ! from toa_radiation.f
c
c ---------------------------------------------------------------------- 
c * * * calendar and orbital calculations * * *
c ---------------------------------------------------------------------- 
c
c calculate time in hours
c
      time = get_time(istep)
      rtime = time / 3600.0
c
c calculate some variables related to top-of-atmosphere radiation
c
      orbit = calc_orbit(jday)
      xdecl = calc_xdecl(orbit)
      sw = calc_sw(orbit)
c
c 9001 continue
c
c do for all gridcells
c
      DO 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c ---------------------------------------------------------------------- 
c * * * solar transmission & diffuse/direct calculations * * *
c ---------------------------------------------------------------------- 
c
c calculate the latitude in radians
c
        jj = latindex(i)
c
        xlat = latscale(jj) * pi / 180.0
c
c fetch the cosine of the solar zenith angle from pre-computed values
c
        coszen(i) = coszen_pre(jj, istep)
c
c find daylength to be used in pheno subroutine
c
c WJS 02.03.10: This equation used to be different, but I believe that
c the old version gave incorrect results at any time step when coszen(i)
c was > 0. This new version gives the same (correct) daylength in each
c time step. Thus it should probably be moved outside of 'diurnal' (and
c outside of the sub-daily loop entirely).
c
        daylength(i) = (180./pi)*((2.*60.)/15.)*
     >    (ACOS(-(SIN(xlat)*SIN(xdecl)) / (COS(xlat)*COS(xdecl))))
c
c calculate the solar transmission through the atmosphere using direct
c measurements
c 
c transmission is fraction of insolation that makes it to the surface
c
        trans(i) = rads(i) / sw
c
c calculate the fraction of indirect (diffuse) solar radiation
c based upon the cloud cover
c
c note that these relationships typically are measured for either
c monthly or daily timescales, and may not be exactly appropriate
c for hourly calculations
c
c method i --
c
c we use a simple empirical relationships from Nikolov and Zeller (1992)
c
c Nikolov, N. and K.F. Zeller, 1992:  A solar radiation algorithm for ecosystem
c dynamics models, Ecological Modelling, 61, 149-168.
c
        fdiffuse = 1.0045 + 0.0435 * trans(i) 
     >                    - 3.5227 * trans(i)**2
     >                    + 2.6313 * trans(i)**3
c
        IF (trans(i).GT.0.75) fdiffuse = 0.166
c
c method ii --
c
c another method was suggested by Spitters et al. (1986) based on
c long-term data from the Netherlands
c
c Spitters et al., 1986: Separating the diffuse and direct component
c of global radiation and its implications for modeling canopy
c photosynthesis, Part I: Components of incoming radiation,
c Agricultural and Forest Meteorology, 38, 217-229.
c
c       if ((trans(i).eq.0.00).and.(trans(i).lt.0.07)) then
c         fdiffuse = 1.0
c       else if ((trans(i).ge.0.07).and.(trans(i).lt.0.35)) then
c         fdiffuse = 1.0 - 2.3 * (trans(i) - 0.07)**2
c       else if ((trans(i).ge.0.35).and.(trans(i).lt.0.75)) then
c         fdiffuse = 1.33 - 1.46 * trans(i)
c       else
c         fdiffuse = 0.23
c       endif
c
c do for each waveband
c
        DO 120 ib = 1, nband
c
c calculate the fraction in each waveband
c
          wfrac = 0.46 + 0.08 * FLOAT(ib - 1)
c
c calculate the direct and indirect solar radiation
c
          solad(i,ib) = wfrac * rads(i) * (1. - fdiffuse) 
c
          solai(i,ib) = wfrac * rads(i) * fdiffuse
c
  120   CONTINUE
c
c ---------------------------------------------------------------------- 
c * * * ir flux calculations * * *
c ---------------------------------------------------------------------- 
c
c clear-sky emissivity as a function of water vapor pressure
c and atmospheric temperature
c
c calculate the ir emissivity of the clear sky
c using equation from Idso (1981) Water Resources Res., 17, 295-304
c
        emb = 0.01 * (psurf(i) * qa(i) / (0.622 + qa(i)))
        ea  = 0.70 + 5.95e-5 * emb * EXP(1500.0 / ta(i))
c
c assign the ir emissivity of clouds (assume they are ~black in the ir)
c
        ec = 0.950
c
c assign the temperature difference of emissions (air + cloud) from
c the surface air temperature
c
        dtair   = 2.0
        dtcloud = 2.0
c
c total downward ir is equal to the sum of:
c
c (1) clear sky contribution to downward ir radiation flux
c (2) cloud contribution to downward ir radiation flux
c
c
        truecloud = 1. - ((trans(i) - 0.251) / 0.509) 
        fira(i) = (1. -  truecloud) * ea * stef * (ta(i) - dtair  )**4 +
     >                   truecloud  * ec * stef * (ta(i) - dtcloud)**4
c
c
c ---------------------------------------------------------------------- 
c * * * snow and rain rate calculations * * *
c ---------------------------------------------------------------------- 
c
c reset snow and rain rate to zero
c
        snowa(i) = 0.0
        raina(i) = 0.0
c
c
c for rain / snow partitioning, make it all rain if 
c ta > 2.5 C and all snow if ta <= 2.5 C
c
c reference:
c
c Auer, A. H., 1974: The rain versus snow threshold temperatures,
c Weatherwise, 27, 67.
c
c
        IF (ta(i)-273.15 .GT. 2.5) THEN
           raina(i) = precip(i) / 3600.0
           IF(i .EQ. 1 .AND. raina(1) .GT. 0) THEN
           ENDIF
        ELSE
           snowa(i) = precip(i) / 3600.0
           IF(i .EQ. 1 .AND. snowa(1) .GT. 0) THEN
            ENDIF
        ENDIF

c
c ---------------------------------------------------------------------- 
c * * * irrigation calculations * * *
c ---------------------------------------------------------------------- 
c
c reset rate of irrigation application per timestep 
c
        xirriga(i) = 0.0
c
c if precipitation event - then no irrigation that day 
c
        IF (time.GE.starti .AND. time.LT.endi
     >      .AND. doirrigate  .GE. 1
     >      .AND. precip(i) .EQ. 0.00) THEN  
c
          xirriga(i) = xirrig(i) / ilens
c
c update annual total - totirrig
c rate of irrigation multiplied by length of timestep (mm/s * s) = total applied
c for this timestep 
c
          totirrig(i) = totirrig(i) + (xirriga(i) * dtime)
c
        ENDIF 
c
  100 CONTINUE
c
      RETURN
      END
c
c -----------------------------------------------------------------------------------
      SUBROUTINE dailyMetCalc (imonth, iday, clear)
c -----------------------------------------------------------------------------------
c
c this routine updates GDDs and such
c
c uses:
c
      USE comgrid
      USE compar
      USE combcs
      USE comatm
      USE comsum
      USE comveg
      USE comsoi
      USE comcrop
c
      IMPLICIT NONE
c
c Arguments
c
      INTEGER imonth,
     >     iday,
c
c local variables
c
     >     i,j,k,               ! loop indice
     >     clear                ! check if gdd needs clearing
c
c
c ---------------------------------------------------------------------- 
c * * * initial setup for daily climate calculations * * *
c ---------------------------------------------------------------------- 
c
c initialize this year's values of gdd0, gdd5, tc, tw at beginning of year
c
      IF (clear .EQ. 1) THEN
c
         CALL CONST (tcthis, npoi,  100.0)
         CALL CONST (twthis, npoi, -100.0)
c     
         CALL CONST (gdd0this, npoi, 0.0)
         CALL CONST (gdd5this, npoi, 0.0)
         CALL CONST (gdd0cthis, npoi, 0.0)
         CALL CONST (gdd8this, npoi, 0.0)
         CALL CONST (gdd10this, npoi, 0.0)
c     
c initialize this year's value of npp and soil/plant gdds (for non multi-year crops)
c
         DO 20 i = 1, npoi
           if(cmask(i) == 0 ) cycle

            DO 50 j = 1, npft
c     
               IF (j .LE. scpft-1) THEN ! natural vegetation  
                  ayanpp(i,j)     = 0.0
               ELSE IF (croplive(i,j) .EQ. 0 .AND. j .GE. scpft) THEN ! non multi-year crops
                  gddplant(i,j)   = 0.0
                  gddtsoi(i,j)    = 0.0
                  ayanpp(i,j)     = 0.0
               ENDIF
c
 50       CONTINUE
 20     CONTINUE
c     
      ENDIF
c
c ---------------------------------------------------------------------- 
c * * * set daily climatic variables for entire domain * * *
c ---------------------------------------------------------------------- 
c
      DO 200 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c 
c calculated temperature extremes -- for vegetation limits (deg c)
c
c for this purpose, use the 10-day running mean temperature
c
         tcthis(i) = MIN (tcthis(i), (a10td(i) - 273.16))
         twthis(i) = MAX (twthis(i), (a10td(i) - 273.16))
c
c update this year's growing degree days
c
         gdd0this(i)  = gdd0this(i)  + MAX (0.0 ,(td(i) - 273.16))   
         gdd0cthis(i) = gdd0cthis(i) + MAX (0.0 ,(td(i) - 273.16)) ! wheat uses base 0C
         gdd5this(i)  = gdd5this(i)  + MAX (0.0 ,(td(i) - 278.16))
         gdd8this(i)  = gdd8this(i)  + MAX (0.0 ,(td(i) - 281.16)) ! maize uses base 8C
         gdd10this(i) = gdd10this(i) + MAX (0.0 ,(td(i) - 283.16)) ! soybean uses base 10C
c     
c accumulate soil/plant growing degree days for planted crops
c
         DO 150 j = scpft, ecpft
c
c for crops except winter wheat
c 
            IF (croplive(i,j) .EQ. 1.0 .AND. iwheat .NE. 2) THEN
c
               gddplant(i,j) = gddplant(i,j) + MAX(0.0, MIN(td(i)
     >              - baset(j), mxtmp(j)))
c
               gddtsoi(i,j)  = gddtsoi(i,j) + MAX(0.0, MIN(tsoi(i,1)
     >              - baset(j), mxtmp(j)))
c
c if winter wheat is planted, reduce thermal time accumulation
c by vernalization factor (calculated in crops.f) until crop
c is fully vernalized (cold-hardened)
c
            ELSE IF (croplive(i,j) .EQ. 1.0 .AND. iwheat .EQ. 2) THEN
c
               gddplant(i,j) = gddplant(i,j) + vf(i) * 
     >              MAX(0.0, MIN(td(i) - baset(j), mxtmp(j)))
c
               gddtsoi(i,j)  = gddtsoi(i,j) +  vf(i) * 
     >              MAX(0.0, MIN(tsoi(i,1) - baset(j), mxtmp(j)))
c
            END IF
c
 150    CONTINUE
c     
 200  CONTINUE
c
c return to main program
c
      RETURN
      END
