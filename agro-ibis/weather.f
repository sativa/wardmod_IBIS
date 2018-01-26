c
c #    #  ######    ##     #####  #    #  ######  #####
c #    #  #        #  #      #    #    #  #       #    #
c #    #  #####   #    #     #    ######  #####   #    #
c # ## #  #       ######     #    #    #  #       #####
c ##  ##  #       #    #     #    #    #  #       #   #
c #    #  ######  #    #     #    #    #  ######  #    #
c
c ---------------------------------------------------------------------
      subroutine daily (iyear, imonth, iday, jday, seed, jdaily, iyrlast, nrun)
c ---------------------------------------------------------------------
c
c overview
c
c this routine generates daily weather conditions from monthly-mean
c climatic parameters
c
c specifically, this routine generates daily values of
c
c  - daily total precipitation
c  - daily maximum temperature
c  - daily minimum temperature
c  - daily average cloud cover
c  - daily average relative humidity
c  - daily average wind speed
c
c in order to generate daily weather conditions, the model uses a series
c of 'weather generator' approaches, which generate random combinations of
c weather conditions based upon the climatological conditions
c
c in general, this weather generator is based upon the so-called Richardson
c weather generator
c
c appropriate references include:
c
c Geng, S., F.W.T. Penning de Vries, and L. Supit, 1985:  A simple
c method for generating rainfall data, Agricultural and Forest
c Meteorology, 36, 363-376.
c
c Richardson, C. W. and Wright, D. A., 1984: WGEN: A model for 
c generating daily weather variables: U. S. Department of
c Agriculture, Agricultural Research Service.
c
c Richardson, C., 1981: Stochastic simulation of daily
c precipitation, temperature, and solar radiation. Water Resources 
c Research 17, 182-190.
c
c uses:
c
      use comgrid
      use compar
      use cominput
      use combcs
      use comatm
      use comsum
      use comveg
      use comsoi
      use comcrop
c
      implicit none
c
c Arguments
c
      integer seed, 
     >        iyear,
     >        imonth,
     >        iday,
     >        jday,
     >        jdaily   ! 1 if reading in daily weather data
     >                 ! 0 if using random/statistical weather generator      
c
c local variables
c
      integer it1w,      ! indice of previous month (interpolation)
     >        it2w,      ! indice of following month (interpolation)
     >        i,j, k     ! loop indice
c
      real rwork,          ! 
     >     omcloud,        ! cloud cover
     >     omqd,           ! humidity
     >     omtmax,         ! maximum temperature
     >     ran2,           ! function random number generator
     >     dt,             ! used for interpolation
     >     pwet,           ! monthly-average probability of rainy day
     >     pwd,            ! probability of a wet day after a dry day
     >     pww,            ! probability of a wet day after a wet day
     >     rndnum,         ! random number to decide if wet or dry day
     >     rainpwd,        ! average rainfall per wet day
     >     alpha,          ! parameter for gamma function
     >     beta,           ! parameter for gamma function
     >     aa,
     >     ab,
     >     tr1,
     >     tr2,
     >     rn1,rn2, rn3,rn,  !random numbers
     >     s1,
     >     s2,
     >     s12,
     >     z,
     >     tdm,            ! mean daily mean temperature
     >     trngm,          ! mean daily mean temperature
     >     tmaxm,          ! mean maximum temperature
     >     tminm,          ! mean minimum temperature
     >     tmaxd,          ! maximum temperatures for dry days
     >     tmaxw,          ! maximum temperatures for wet days
     >     tmaxe,          !'expected' maximum temperature for today
     >     tmine,          !'expected' minimum temperature for today
     >     tmaxs,          ! standard deviation in maximum temperature (K)
     >     tmins,          ! standard deviation in minimum temperature (K)
     >     cloudm,         ! mean cloud cover for today (fraction)
     >     cloudd,         ! dry day cloud cover
     >     cloudw,         ! wet day cloud cover
     >     cloude,         ! expected cloud cover today
     >     clouds,         ! standard deviation of cloud fraction
     >     v,
     >     tdum,           ! storage variable
     >     qdm,            ! mean relative humidity
     >     qdd,            ! dry day relative humidity
     >     qdw,            ! wet day relative humidity 
     >     qde,            ! expected relative humidity (based on wet/dry decision)
     >     qdup,           ! upper bound of humidity distribution function
     >     qdlow,          ! lower bound of humidity distribution function
     >     y,
     >     b3,
     >     b2,
     >     b1,
     >     x1,
     >     amn,
     >     eud             ! expected daily average wind speed from monthly mean
 
c
      real a(3,3),
     >     b(3,3)
c
      real
     >     ee(3), 
     >     r(3),
     >     rr(3),
     >     x(3)
c
      integer
     >     iyrlast,
     >     nrun
c
      real
     >     precipfac,
     >     dif
c
c Externals
c
      real gdd_natveg, gdd_crop, gdd_soil  ! from gdd.f
c
c define autocorrelation matrices for Richardson generator
c
c note that this matrix should be based upon a statistical
c analysis of regional weather patterns
c
c for global simulations, we use 'nominal' values
c
      data a /  0.600,  0.500,  0.005,
     >          0.010,  0.250,  0.005, 
     >          0.020,  0.125,  0.250 /
c
      data b /  0.500,  0.250, -0.250,
     >          0.000,  0.500,  0.250, 
     >          0.000,  0.000,  0.500 /
c
      include 'comsat.h'
c
c ---------------------------------------------------------------------- 
c * * * initial setup for daily climate calculations * * *
c ---------------------------------------------------------------------- 
c
c define working variables
c
      rwork = (grav / rair / 0.0065)
c
c 'omega' parameters used to calculate differences in expected
c climatic parameters on wet and dry days
c
c following logic of weather generator used in the EPIC crop model
c
c omcloud -- cloud cover
c omqd    -- humidity
c omtmax  -- maximum temperature
c
      omcloud = 0.90    ! originally 0.90
      omqd    = 0.50    ! originally 0.50
      omtmax  = 0.75    ! originally 0.75
c
c calculate weighting factors used in interpolating climatological
c monthly-mean input values to daily-mean values
c
c this is a simple linear interpolation technique that takes into
c account the length of each month
c
      if (jdaily .eq. 0) then
        if (float(iday).lt.float(ndaypm(imonth)+1)/2.0) then
          it1w = imonth - 1
          it2w = imonth 
          dt   = (float(iday) - 0.5) / ndaypm(imonth) + 0.5
        else
          it1w = imonth 
          it2w = imonth + 1
          dt   = (float(iday) - 0.5) / ndaypm(imonth) - 0.5
        end if
c
        if (it1w.lt. 1) it1w = 12
        if (it2w.gt.12) it2w = 1
      end if
c
c initialize this year's values of gdd0, gdd5, tc, tw
c
      if (iday.eq.1 .and. imonth.eq.1) then
c
        call const (tcthis, npoi,  100.0)
        call const (twthis, npoi, -100.0)
c
        call const (gdd0this, npoi, 0.0)
        call const (gdd5this, npoi, 0.0)
        call const (gdd0cthis, npoi, 0.0)
        call const (gdd8this, npoi, 0.0)
        call const (gdd10this, npoi, 0.0)
c
        call const (consdays, npoi, 0.)
        call const (iniday, npoi, 9999.)
        call const (maxcons, npoi, 0.)
        call const (gsdays, npoi, 0.)
        call const (gddfzcorn, npoi, 0.0)
        call const (gddfzsoy, npoi, 0.0)
c    
c
c
c initialize variables to zero at beginning of year
c for crop types that do not grow over two calendar years
c 
c constant for gdd based on 10 cm soil temperature (planting bed-
c used to calculate leaf emergence after planting
c
        do 20 i = 1, npoi
          if(cmask(i) == 0 ) cycle

          do 50 j = 1, npft
c
            if (j .le. scpft-1) then  ! natural vegetation 
              ayanpp(i,j)     = 0.0
            else if (croppresent(i,j) .eq. 0 .and. j .ge. scpft) then
              gddplant(i,j)   = 0.0
              gddtsoi(i,j)    = 0.0
              ayanpp(i,j)     = 0.0
            endif
 50       continue
 20     continue
c           
      endif

c ---------------------------------------------------------------------- 
c make sure we're not trying to use the weather generator with incompatible options
c ---------------------------------------------------------------------- 

      if (jdaily .eq. 0) then

        if (read_radiation) then
          write(*,*) 'ERROR IN WEATHER.F: DAILY:'
          write(*,*) 'Weather generator not implemented for read_radiation = TRUE'
          stop 1
        end if

c Check if the first xinwet point in this month is equal to un_init
c (i.e., uninitialized), within rounding error. If so, this signals that
c we didn't read in any wetd data (which are only needed for the weather
c generator).
        if (abs(xinwet(1, imonth) - un_init) .lt. (epsilon * abs(un_init))) then
          write(*,*) 'ERROR IN WEATHER.F: DAILY:'
          write(*,*) 'No wetd data read in; these are needed for the weather generator'
          stop 1
        end if

      end if  ! if (jdaily .eq. 0)

c ---------------------------------------------------------------------- 
c * * * use weather generator to create daily statistics * * *
c ---------------------------------------------------------------------- 
c
      if (jdaily .eq. 0) then

        do 200 i = 1, npoi
          if(cmask(i) == 0 ) cycle

c ---------------------------------------------------------------------- 
c (1) determine if today will rain or not (following Geng et al.)
c ---------------------------------------------------------------------- 
c
c implement simple first-order Markov-chain precipitation generator logic
c based on Geng et al. (1986), Richardson and Wright (1984),
c and Richardson (1981) 
c
c basically, this allows for the probability of today being a wet day 
c (a day with measureable precipitation) to be a function of what
c yesterday was (wet or dry)
c
c the logic here is that it is more likely that a wet day will follow
c another wet day -- allowing for 'storm events' to persist
c
c calculate monthly-average probability of rainy day 
c
          pwet = max (1., xinwet(i,imonth)) / ndaypm(imonth)
c
c estimate the probability of a wet day after a dry day
c
          pwd = 0.75 * pwet
c
c estimate the probability of a wet day after a wet day
c
          pww = 0.25 + pwd
c
c Beginning of block of code that figures out daily precip for
c entire month on the first day of each month
c
          if (iday .eq. 1) then
c
c Verify the dataset consistency especially when using interannual anomalies
c of precipitations (negative values, or too few wet days in a rainy month)
c
            xinprec(i, imonth) = max(0., xinprec(i, imonth))
            xinwet(i, imonth) = max(1., xinwet(i, imonth))
c
 9000       continue
c
c Initialize monthly parameters back to zero
c
            iwetdaysum(i) = 0
            precipdaysum(i) = 0
c
            do 210 j = 1, 31
              iwetday(i,j) = 0
              precipday(i,j) = 0
 210        continue
c
c Loop through number of days in this month and determine precip
c
            do 220 j = 1, ndaypm(imonth)
c
c decide if today is a wet day or a dry day using a random number
c
              rndnum = ran2(seed)
c
c
c If it is the first day of the month do not look at previous day
c
              if (j .eq. 1) then
                if (rndnum .le. pwd) then
                  iwetday(i,j) = 1
                  iwetdaysum(i) = iwetdaysum(i) + 1
                else
                  iwetday(i,j) = 0
                endif
c
              else
c
c If it is not the first day, look at yesterday's wet/dry index to help
c determine if today is wet or dry
c
                if (iwetday(i,j-1) .eq. 0) then
                  if (rndnum.le.pwd) then
                    iwetday(i,j) = 1
                    iwetdaysum(i) = iwetdaysum(i) + 1
                  endif
                else
                  if (rndnum.gt.pww) iwetday(i,j) = 0
                endif
c
              endif
c
c ---------------------------------------------------------------------- 
c (2) determine today's precipitation amount (following Geng et al.)
c ---------------------------------------------------------------------- 
c
c if it is going to rain today
c
              if (iwetday(i,j) .eq. 1) then
c
c calculate average rainfall per wet day
c
                rainpwd = xinprec(i,imonth) * ndaypm(imonth) /
     >                max (0.1, xinwet(i,imonth))
c
c randomly select a daily rainfall amount from a probability density
c function of rainfall
c
c method i --
c
c use the following technique from Geng et al. and Richardson
c to distribute rainfall probabilities
c 
c pick a random rainfall amount from a two-parameter gamma function
c distribution function
c
c estimate two parameters for gamma function (following Geng et al.)
c 
                beta  = max (1.0, -2.16 + 1.83 * rainpwd)
                alpha = rainpwd / beta
c
c determine daily precipitation amount from gamma distribution function
c (following WGEN code of Richardson and Wright (1984))
c
                aa = 1.0 / alpha 
                ab = 1.0 / (1.0 - alpha)
c
                tr1 = exp(-18.42 / aa)
                tr2 = exp(-18.42 / ab)
c
 12             rn1 = ran2(seed)
                rn2 = ran2(seed)
c
c CD: rewrote parts of prehistoric code in fortran 77 
c
                if ((rn1 - tr1) .le. 0) then
                  s1 = 0.0
                else 
                  s1 = rn1**aa
                end if
c
                if ((rn2 - tr2) .le. 0) then 
                  s2 = 0.0
                else 
                  s2 = rn2**ab
                end if
           
c               if (rn1 - tr1) 61, 61, 62
c 61            s1 = 0.0
c               go to 63
c 62            s1 = rn1**aa
c
c 63            if (rn2 - tr2) 64, 64, 65
c 64            s2 = 0.0
c               go to 66
c 65            s2 = rn2**ab
c
c
 66             s12 = s1 + s2
c
                if (s12 - 1.0)  13, 13, 12
 13             z = s1 / s12
c
                rn3 = ran2(seed)
c
                precipday(i,j) = -z * log(rn3) * beta
c
c method ii --
c
c here we use a one-parameter Weibull distribution function
c following the analysis of Selker and Haith (1990)
c
c Selker, J.S. and D.A. Haith, 1990: Development and testing of single-
c parameter precipitation distributions, Water Resources Research, 
c 11, 2733-2740.
c
c this technique seems to have a significant advantage over other
c means of generating rainfall distribution functions
c
c by calibrating the Weibull function to U.S. precipitation records,
c Selker and Haith were able to establish the following relationship
c
c the cumulative probability of rainfall intensity x is given as:
c
c f(x) = 1.0 - exp(-(1.191 x / rainpwd)**0.75)
c 
c where x       : rainfall intensity
c       rainpwd : rainfall per wet day
c
c using transformation method, take uniform deviate and convert it to a
c random number weighted by the following Weibull function
c
c          rndnum = ran2(seed)
c
c          precip(i) = rainpwd / 1.191 * (-log(1.0 - rndnum))**1.333333
c 
c bound daily precipitation to "realistic" range
c 
c lower end is determined by definition of a 'wet day' (at least
c 0.25 mm of total precipitation)
c
c upper end is to prevent ibis from blowing up
c
                precipday(i,j) = max (precipday(i,j),0.25) ! min =   0.25 mm/day
                precipday(i,j) = min (precipday(i,j),150.00) ! max = 150.00 mm/day
c
c Back to beginning of month loop, this is the end of it
c
              endif
c
c Add today's precip to the monthly summation
c
c
              precipdaysum(i) = precipdaysum(i) + precipday(i,j)
c
 220        continue
c
c Adjust daily precip amounts (using precipfac) so that the monthly
c summation equals the input precip amount, when using interannual
c anomalies
c
c
              if ((precipdaysum(i) .eq. 0) .AND.
     >                (xinprec(i,imonth) .gt. 0.)) then
                    rndnum = 1.0 + (float(ndaypm(imonth)) - 1.0)
     >                   * ran2(seed)
                    iwetday(i,nint(rndnum)) = 1
                    precipday(i,nint(rndnum)) = xinprec(i,imonth)
     >                   * float(ndaypm(imonth))
                    precipdaysum(i) = precipday(i,nint(rndnum))
                    iwetdaysum(i) = 1
                 end if
c
                 precipfac = (xinprec(i,imonth)*float(ndaypm(imonth)))
     >                / max(0.01,precipdaysum(i))
c
                 do 230 j=1,ndaypm(imonth)
                    precipday(i,j) = precipday(i,j) * precipfac
c
                    if (precipday(i,j).gt.360) then
                       if (xinwet(i,imonth) .lt. ndaypm(imonth)) then
                          xinwet(i,imonth) = xinwet(i,imonth) + 1
                          pwet = xinwet(i,imonth) / ndaypm(imonth)
                          pwd = 0.75 * pwet
                          pww = 0.25 + pwd
                          print *, 'WARNING: goto 9000a', i, int(xinwet(i
     >                         ,imonth)), iwetdaysum(i), int(precipday(i,j))
                          goto 9000
                       else
                          print *, 'WARNING: goto 9000b', i, int(xinwet(i
     >                         ,imonth)), iwetdaysum(i), int(precipday(i,j))
                          goto 9000
                       end if
                    end if
c    
  230            continue
c    
c Verification of the weather generator algorithm
c
                 iwetdaysum(i) = 0
                 precipdaysum(i) = 0.
c
                 do 240 j=1,ndaypm(imonth)
                    precipdaysum(i) = precipdaysum(i) + precipday(i,j)
                    iwetdaysum(i) = iwetdaysum(i) + iwetday(i,j)
  240            continue
c
                 dif = precipdaysum(i) - xinprec(i,imonth)
     >                   * float(ndaypm(imonth))
c    
                 if ((dif.lt.-0.1).or.(dif.gt.0.1)) print *,
     >                'ERROR in DAILY:', i, precipdaysum(i),
     >                xinprec(i,imonth)* float(ndaypm(imonth)),
     >                iwetdaysum(i), xinwet(i,imonth)
c    
c end of the verification
c
           end if               !end of the iday loop
c    
c Relate today's iwetday and precipday to iwet and precip that will be
c used below
           iwet(i) = iwetday(i,iday)
           precip(i) = precipday(i,iday)
c
c ---------------------------------------------------------------------- 
c (3) estimate expected minimum and maximum temperatures
c ---------------------------------------------------------------------- 
c
c first determine the expected maximum and minimum temperatures
c (from climatological means) for this day of the year
c
c mean daily mean temperature (K)
c
          tdm = xint(i,it1w) +
     >          dt * (xint(i,it2w) - xint(i,it1w)) + 273.16
c
c mean daily temperature range (K)
c
          trngm = xintrng(i,it1w) +
     >            dt * (xintrng(i,it2w) - xintrng(i,it1w))
c
c mean minimum and maximum temperatures
c
          tmaxm = tdm + tmin_wt * trngm
          tminm = tdm - tmax_wt * trngm
c
c modify maximum temperatures for wet and dry days
c
          if (pwet .ne. 0.0) then
            tmaxd = tmaxm + pwet * omtmax * trngm
            tmaxw = tmaxd -        omtmax * trngm
          else
            tmaxd = tmaxm
            tmaxw = tmaxm
          endif
c
c set the 'expected' maximum and minimum temperatures for today
c
c note that the expected minimum temperatures are the same for
c both wet and dry days
c
          if (iwet(i).eq.0) tmaxe = tmaxd
          if (iwet(i).eq.1) tmaxe = tmaxw
c
          tmine = tminm
c
c estimate variability in minimum and maximum temperatures
c
c tmaxs : standard deviation in maximum temperature (K)
c tmins : standard deviation in minimum temperature (K)
c
c Regression is based on analysis of 2-m air temperature data from the
c NCEP/NCAR reanalysis (1958-1997) for 294 land points over central
c North America (24N-52N, 130W-60W, 0.5-degree resolution): Daily max
c and min temperatures were calculated for each land point from daily
c mean temperature and temperature range. Anomalies were calculated
c by subtracting similar max and min temperatures calculated from
c monthly mean temperature and range (interpolated to daily). The 40
c years of anomalies were then binned by month and the standard
c deviation calculated for each month. The 294 x 12 standard
c deviations were then regressed against the 3528 long-term monthly
c mean temperatures.
c
c note: the values are bound to be greater than 1.0 K 
c       (at the very least they must be bound so they don't go below zero)
c
          tmaxs = max (1.0, -0.0713 * (tdm - 273.16) + 4.89)
          tmins = max (1.0, -0.1280 * (tdm - 273.16) + 5.73)
c
c ---------------------------------------------------------------------- 
c (4) estimate expected cloud cover
c ---------------------------------------------------------------------- 
c
c the formulation of dry and wet cloud cover has been
c derived from the weather generator used in the epic crop model
c
c cloudm : mean cloud cover for today
c cloudd : dry day cloud cover
c cloudw : wet day cloud cover
c cloude : expected cloud cover today
c
c Verify the data set consistency when using interannual anomalies of
c cloudiness (values under 0 % or over 100 %)
c
          if (iday.eq.1) then
             xincld(i,it1w) = max (0.0, xincld(i,it1w))
             xincld(i,it1w) = min (100., xincld(i,it1w))
             xincld(i,it2w) = max (0.0, xincld(i,it2w))
             xincld(i,it2w) = min (100., xincld(i,it2w))
          end if

c monthly mean cloud cover (%)
c
          cloudm = xincld(i,it1w) +
     >             dt * (xincld(i,it2w) - xincld(i,it1w))
c 
c convert from percent to fraction
c
          cloudm = cloudm / 100.0
c
c adjust cloud cover depending on dry day / rainy day
c following logic of the EPIC weather generator code
c
          if (pwet .ne. 0.0) then
            cloudd = (cloudm - pwet * omcloud) / (1.0 - pwet * omcloud)
            cloudd = min (1.0, max (0.0, cloudd))
            cloudw = (cloudm - (1.0 - pwet) * cloudd) / pwet
          else
            cloudd = cloudm
            cloudw = cloudm
          endif
c
          if (iwet(i).eq.0) cloude = cloudd
          if (iwet(i).eq.1) cloude = cloudw
c
c estimate variability in cloud cover for wet and dry days
c following numbers proposed by Richardson
c
c clouds : standard deviation of cloud fraction
c 
          if (iwet(i).eq.0) clouds = 0.24 * cloude
          if (iwet(i).eq.1) clouds = 0.48 * cloude
c
c ---------------------------------------------------------------------- 
c (5) determine today's temperatures and cloud cover using
c     first-order serial autoregressive technique
c ---------------------------------------------------------------------- 
c
c use the Richardson (1981) weather generator approach to simulate the
c daily values of minimum / maximum temperature and cloud cover
c
c following the implementation of the Richardson WGEN weather generator
c used in the EPIC crop model
c
c this approach uses a multivariate generator, which assumes that the
c perturbation of minimum / maximum temperature and cloud cover are
c normally distributed and that the serial correlation of each
c variable may be described by a first-order autoregressive model
c
c generate standard deviates for weather generator
c
          do j = 1, 3
 31         rn1 = ran2(seed)
            rn2 = ran2(seed)
            v = sqrt (-2.0 * log(rn1)) * cos(6.283185 * rn2)
            if (abs(v) .gt. 2.5) go to 31
            ee(j) = v
          enddo
c
c zero out vectors
c
          do j = 1, 3
            r(j)  = 0.0
            rr(j) = 0.0
          enddo
c
c update working vectors
c
          do j = 1, 3
            do k = 1, 3
              r(j)  = r(j)  + b(j,k) * ee(j)
              rr(j) = rr(j) + a(j,k) * xstore(i,k)
            enddo
          enddo
c
c solve for x() perturbation vector and save current vector
c into the xim1() storage vector (saved for each point)
c 
          do j = 1, 3
            x(j) = r(j) + rr(j)
            xstore(i,j) = x(j)
          enddo
c
c determine today's minimum and maximum temperature
c
          tmax(i)  = tmaxe + tmaxs * x(1)
          tmin(i)  = tmine + tmins * x(2)
c
c if tmin > tmax, then switch the two around
c
          if (tmin(i).gt.tmax(i)) then
            tdum    = tmax(i)
            tmax(i) = tmin(i)
            tmin(i) = tdum
          endif
c
c daily average temperature
c
          td(i) = tmax_wt * tmax(i) + tmin_wt * tmin(i)
c
c determine today's cloud cover
c
          cloud(i) = cloude + clouds * x(3)
c
c constrain cloud cover to be between 0 and 100%
c
          cloud(i) = max (0.0, min (1.0, cloud(i)))
c
c estimate transmission fraction and solar radiation at the ground based
c on cloud cover
c
          call derived_vars_from_cloud (jday)
c
c ---------------------------------------------------------------------- 
c (6) estimate today's surface atmospheric pressure
c ---------------------------------------------------------------------- 
c
c simply a function of the daily average temperature and topographic
c height -- nothing fancy here
c
          psurf(i) = 101325.0 *
     >               (td(i) / (td(i) + 0.0065 * xintopo(i))) ** rwork
c
c ---------------------------------------------------------------------- 
c (7) estimate today's relative humidity
c ---------------------------------------------------------------------- 
c
c the formulation of dry and wet relative humidities has been
c derived from the weather generator used in the epic crop model
c
c qdm : mean relative humidity 
c qdd : dry day relative humidity
c qdw : rainy day relative humidity
c qde : expected relative humidity (based on wet/dry decision)
c
c Verify the data set consistency when using interannual anomalies of
c relative humidity (values over 100 % or under 0 %)
c
          if (iday.eq.1) then
             xinq(i,it1w) = max (0.0, xinq(i,it1w))
             xinq(i,it1w) = min (100., xinq(i,it1w))
             xinq(i,it2w) = max (0.0, xinq(i,it2w))
             xinq(i,it2w) = min (100., xinq(i,it2w))
          end if
c
c mean relative humidity (%)
c
          qdm = xinq(i,it1w) + dt * (xinq(i,it2w) - xinq(i,it1w))
c 
c convert from percent to fraction
c
          qdm = qdm / 100.0
c
c adjust humidity depending on dry day / rainy day
c following logic of the EPIC weather generator code
c
          if (pwet .ne. 0.0) then
            qdd = (qdm - pwet * omqd) / (1.0 - pwet * omqd)
            if (qdd .lt. 0.2) then
              qdd = 0.2
              if (qdd .gt. qdm) qdm = qdd
            endif
            qdd = min(1.0, qdd)
            qdw = (qdm - (1.0 - pwet) * qdd) / pwet
          else
            qdd = qdm
            qdw = qdm
          endif
c
          if (iwet(i).eq.0) qde = qdd
          if (iwet(i).eq.1) qde = qdw 
c
c estimate lower and upper bounds of humidity distribution function
c following logic of the EPIC weather generator code
c
          qdup  = qde + (1.0 - qde) * exp (qde - 1.0)
          qdlow = qde * (1.0 - exp (-qde))
c
c randomly select humidity from triangular distribution function
c following logic of the EPIC weather generator code
c
          rn = ran2(seed)
c
          y  = 2.0 / (qdup - qdlow)
c
          b3 = qde  - qdlow
          b2 = qdup - qde
          b1 = rn / y
c
          x1 = y * b3 / 2.0 
c
          if (rn.gt.x1) then
            qd(i) = qdup  - sqrt (b2 * b2 - 2.0 * b2 * (b1 - 0.5 * b3))
          else
            qd(i) = qdlow + sqrt (2.0 * b1 * b3)
          endif
c
c adjust daily humidity to conserve monthly mean values
c
c note that this adjustment sometimes gives rise to humidity
c values greater than 1.0 -- which is corrected below
c
          amn = (qdup + qde + qdlow) / 3.0
          qd(i) = qd(i) * qde / amn
c
c constrain daily average relative humidity
c
          qd(i) = max (0.30, qd(i))
          qd(i) = min (0.99, qd(i))
c
c convert from relative humidity to specific humidity at
c daily mean temperature
c
          qd(i) = qsat(qd(i) * esat(td(i)), psurf(i))
c
c ---------------------------------------------------------------------- 
c (8) estimate today's daily average wind speed
c ---------------------------------------------------------------------- 
c
c first estimate the expected daily average wind speed (from monthly
c means)
c
          eud = xinwind(i,it1w) +
     >          dt * (xinwind(i,it2w) - xinwind(i,it1w))
c
c following logic of the EPIC weather generator
c select random wind speed following this equation
c
          ud(i) = 1.13989 * eud * (-log(ran2(seed)))**0.30 
c
c WJS (02.10.10): Commented out the following adjustment because:
c  (1) It gave infinite wind speeds for crops (which have ztop(i,2)==0)
c  (2) There was supposed to be a later adjustment back down to a lower
c      height, but it looks like this readjustment never happened.
c So I got rid of this adjustment and instead changed siga back to
c 0.999, which is the approximate height of the wind speed data.
c
c$$$c TET change 10/29/02
c$$$c Input windspeed is assumed to be measured at canopy height
c$$$c (ztop(i,2))
c$$$c IBIS now has siga = 0.991 which is ~93m height.
c$$$c Adjust windspeed according to log-wind profile.
c$$$c So let
c$$$c displacement height, d = 0.65(ztop(i,2))
c$$$c roughness length, zm = 0.1(ztop(i,2))
c$$$c Equation found in Bonan (2002) and other texts:
c$$$c u(z2)-u(z1)=(frictionvelocity/k)*log((z2-d)/(z1-d))
c$$$c Use log-wind to solve for frictionvelocity and substitute
c$$$c into above equation:
c$$$c u(z2) = u(z1)*(1.+(log((z2-d)/(z1-d))) / (log((z1-d)/zm)))
c$$$c Use canopyheight = z1 = ztop(i,2), and z2=93m
c$$$c and substitute d and zm to get equation dependent on canopy height:
c$$$c
c$$$c WJS (01.21.10): Note that 1/(log((z1-d)/zm)) = 0.79824, regardless of
c$$$c the value of ztop(i,2); also, (z1 - d) = 0.35*ztop(i,2)
c$$$
c$$$          ud(i) = ud(i)*(1. + 0.79824*log((93.-0.65*ztop(i,2))/
c$$$     *                   (0.35*ztop(i,2))))
c$$$c
c$$$c Decided to use a canopy height of 20m over US after finding
c$$$c references of mature forest height between 10m and 35m (we are
c$$$c trying to simulate forests before human interference)
c$$$c with ztop(i,2) = 20m, the input windspeed is adjusted by a
c$$$c factor of 3.21
c
c constrain daily wind speeds to be between 2.5 and 25.0 m/sec
c
c WJS 02.10.10: The original constraint was 2.5 - 10, then TET changed
c it to 2.5 - 50 (presumably because ud now referred to wind speed at
c 93m); based on discussions with CJK, we think 10 is too low but 50 is
c too high (now that we no longer do TET's correction to 93 m), so we're
c using a max of 25 (somewhat arbitrarily).
c
c Original:
c          ud(i) = max (2.5, min (10.0, ud(i)))
c TET changed constraint since it will be above 10 often:
c          ud(i) = max (2.5, min (50.0, ud(i)))
c
          ud(i) = max (2.5, min (25.0, ud(i)))

 200    continue

c ---------------------------------------------------------------------- 
c * * * use real daily climate data * * *
c ---------------------------------------------------------------------- 
c
      else

        call create_daily_precips (imonth)
        call create_daily_temps (imonth)
        call create_daily_cloud_and_rads (imonth, jday)

c compute surface atmospheric pressure
        do i = 1, npoi
          if(cmask(i) == 0 ) cycle

          psurf(i) = 101325.0 *
     >               (td(i) / (td(i) + 0.0065 * xintopo(i))) ** rwork
        end do

        call create_daily_humidity (imonth)
        call create_daily_wind (imonth)

      end if
c
c ---------------------------------------------------------------------- 
c * * * other daily climate calculations * * *
c ---------------------------------------------------------------------- 
c

      do 202 i = 1, npoi
        if(cmask(i) == 0 ) cycle

c calculated temperature extremes -- for vegetation limits (deg c)
c 
c for this purpose, use the 10-day running mean temperature
c 
        tcthis(i) = min (tcthis(i), (a10td(i) - 273.16))
        twthis(i) = max (twthis(i), (a10td(i) - 273.16))
c 
c update this year's growing degree days
c 
        gdd0this(i)  = gdd0this(i)  + gdd_natveg(td(i)-c_to_k, 0.0)
        gdd0cthis(i) = gdd0cthis(i) + gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, 0.0, 26.0)  ! wheat
        gdd5this(i)  = gdd5this(i)  + gdd_natveg(td(i)-c_to_k, 5.0)
        gdd8this(i)  = gdd8this(i)  + gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, 8.0, 30.0)  ! maize (prior to 4-22-10)
        gdd10this(i) = gdd10this(i) + gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, 10.0, 30.0) ! soybean (& maize, since 4-22-10)
c 
c form calculations of number of growing degree days between frost events
c 
c events (e.g., tmin(i) .le. -2.2 C) this is applicable to CRU data
c differences exist when using a combination of CRU/NCEP  
c -2.2 used as a threshold from Schwartz and Reiter, 2000.  International
c Journal of Climatology, 20: 929-932.  Changes in North American Spring 
c 
        if (tmin(i) .ge. 273.16) then
          consdays(i) = consdays(i) + 1
          maxcons(i)  = max(maxcons(i), consdays(i))
          if (maxcons(i) .eq. consdays(i)) then
            iniday(i) = jday - maxcons(i)
          endif
          daygddc(i,jday) = gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(maize)-c_to_k, mxtmp(maize))
          daygdds(i,jday) = gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(soy)-c_to_k, mxtmp(soy))
          daygddw(i,jday) = gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(wheat)-c_to_k, mxtmp(wheat))
c 
        else
          consdays(i) = 0
c 
        endif
c 
        if (iday .eq. 31 .and. imonth .eq. 12) then
c 
c calculate total number of growing season GDD and number
c of days long
c 
          if (iniday(i) .eq. 9999) then
            iniday(i) = 1
            maxcons(i) = jday-1
          endif
          
          endday(i) = iniday(i) + maxcons(i)-1
          do 125 k = iniday(i), iniday(i)+maxcons(i)-1
            gddfzcorn(i) =  gddfzcorn(i) + daygddc(i,k)
            gddfzsoy(i)  =  gddfzsoy(i)  + daygdds(i,k)
            gddfzwht(i)  =  gddfzwht(i)  + daygddw(i,k)
            gsdays(i) = gsdays(i) + 1 
 125      continue
          gddcorn(i,iyear) = gddfzcorn(i) 
          gddsoy(i,iyear)  = gddfzsoy(i) 
          gddwht(i,iyear)  = gddfzwht(i) 
        endif
c 
c accumulate growing degree days for planted crops past planting
c 
        do 150 j = scpft, ecpft 
c 
c for crops except winter wheat
c be careful with rotations
c 
c We use croppresent rather than croplive, so we continue to accumulate GDD between maturity & harvest,
c in case any logic depends on that continued accumulation of GDD
c
          if (croppresent(i,j) .eq. 1 .and. ((j .eq. 13 .or. j .eq. 14)
     >      .or.    (iwheat .eq. 1 .and. j .eq. 15))) then
c 
            gddplant(i,j) = gddplant(i,j) + 
     >        gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c 
            gddtsoi(i,j)  = gddtsoi(i,j) + 
     >        gdd_soil(tsoi(i,1)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c 
c if winter wheat is planted, reduce thermal time accumulation
c by vernalization factor (calculated in crops.f) until crop
c is fully vernalized (cold-hardened)
c 
c We use croppresent rather than croplive, so we continue to accumulate GDD between maturity & harvest,
c in case any logic depends on that continued accumulation of GDD
c
          else if (croppresent(i,j) .eq. 1 .and. j .eq. 15 .and. 
     >        iwheat .eq. 2) then
c 
            gddplant(i,j) = gddplant(i,j) + 
     >        vf(i) * gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c 
            gddtsoi(i,j)  = gddtsoi(i,j)  +  
     >        vf(i) * gdd_soil(tsoi(i,1)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c 
          endif
c 
 150    continue
c 
 202  continue
c
c return to main program
c
      return
      end

c ---------------------------------------------------------------------
c HELPER SUBROUTINES CALLED BY SUBROUTINE 'DAILY'
c ---------------------------------------------------------------------

c The following subroutines set the values of various climate variables
c when we're using daily data (not the weather generator). A single call
c to the subroutine sets the value at all points (for i = 1 to npoi).
c
c For each variable that's read in, there are two options: If
c read_daily_directly is true, then we simply set the variable in the
c model to the value that was read in (e.g., precip(i) = xinprecd(i)),
c possibly converting units. If read_daily_directly is false, then we
c calculate the variable in the model using monthly averages and daily
c anomalies. Note that xin*d stands for the variable that's read in from
c the daily files; this can be either the daily value itself (if
c read_daily_directly is true) or the daily anomalies (if
c read_daily_directly is false).

c ---------------------------------------------------------------------
      subroutine create_daily_precips (imonth)
c ---------------------------------------------------------------------
c
c This subroutine sets precip(i) (for i = 1 to npoi) when using real
c daily data
c
      use comgrid
      use compar
      use cominput
      use combcs
      use comatm
c
      implicit none
c
c Arguments
c
      integer imonth
c
c Local variables
c
      integer i
c
      do i = 1, npoi
        if(cmask(i) == 0 ) cycle

        if (read_daily_directly) then
          precip(i) = xinprecd(i)
        else
c Here we multiply xinprecd, the daily fraction of precip calculated from
c the NCEP dataset, by xinprec, the total monthly amount of precip taken from
c the CRU05 dataset to obtain our derived daily precip amount. Also do a check
c to see if the daily precip exceeds 360mm (as is done in the daily weather 
c generator) ... no correction is made, only a warning is printed
c 
          precip(i) = (xinprec(i,imonth) * ndaypm(imonth)) * xinprecd(i)
c 
c if (precip(i) .gt. 360) then
c print *, 'WARNING: daily precip exceeds 360mm for'
c print *, 'year, month, day, gridcell = '
c print *, iyear, imonth, iday, i
c endif
        end if

      end do

      return
      end
      

c ---------------------------------------------------------------------
      subroutine create_daily_temps (imonth)
c ---------------------------------------------------------------------
c
c This subroutine sets td(i), tmax(i) and tmin(i) (for i = 1 to npoi)
c when using real daily data
c
      use comgrid
      use compar
      use cominput
      use combcs
      use comatm
c
      implicit none
c
c Arguments
c
      integer imonth
c
c Local variables
c
      integer i
      real trngm  ! daily temperature range
c
      if (read_tmin_tmax) then
c Note: For this option, we don't impose a max of 44 on the daily temperature range, as we do in the 'else' clause
        do i = 1, npoi
          if(cmask(i) == 0 ) cycle

          if (read_daily_directly) then
            tmax(i) = xintmaxd(i) + 273.16
            tmin(i) = xintmind(i) + 273.16
          else
            tmax(i) = xintmax(i,imonth) + 273.16 + xintmaxd(i)
            tmin(i) = xintmin(i,imonth) + 273.16 + xintmind(i)
          end if

          td(i) = tmax_wt * tmax(i) + tmin_wt * tmin(i)

        end do

      else  ! .not. read_tmin_tmax
        do i = 1, npoi
          if(cmask(i) == 0 ) cycle

          if (read_daily_directly) then
            td(i) = xintd(i) + 273.16
            trngm = min(44.0, xintrngd(i))
          else
c Here we add the NCEP temperature anomaly to the CRU05 monthly anomaly
c The trange NCEP anomaly must also be multiplied by the climatological
c CRU05 trange in order to determine tmax and tmin
            td(i) = xint(i,imonth) + 273.16 + xintd(i)
            trngm = min (44.0, (xintrng(i,imonth) * xintrngd(i)))
          end if
c 
          tmax(i) = td(i) + tmin_wt * trngm
          tmin(i) = td(i) - tmax_wt * trngm

        end do
      end if  ! read_tmin_tmax
c 
      return
      end


c ---------------------------------------------------------------------
      subroutine create_daily_cloud_and_rads (imonth, jday)
c ---------------------------------------------------------------------
c
c This subroutine sets cloud(i), rads(i) and trans(i) (for i = 1 to
c npoi) when using real daily data
c
      use comgrid
      use compar
      use cominput
      use combcs
      use comatm
c
      implicit none
c
c Arguments
c
      integer imonth,
     >        jday     ! day of the year
c
c Local variables
c
      integer i

      if (read_radiation) then
c Use read-in values of radiation to set rads(i); compute trans(i) and cloud(i) from this        
        do i = 1, npoi
          if(cmask(i) == 0 ) cycle

          if (read_daily_directly) then
            rads(i) = xinradsd(i)
          else
            write(*,*) "ERROR: IN CREATE_DAILY_CLOUD_AND_RADS:"
            write(*,*) "The combination of read_radiation = TRUE"
            write(*,*) "and read_daily_directly = FALSE"
            write(*,*) "is currently unimplemented"
            stop 1
          end if
        end do

        call derived_vars_from_rads (jday) 
      
      else  ! read_radiation = FALSE
c Use read-in values of cloudiness to set cloud(i); compute trans(i) and rads(i) from this        
        do i = 1, npoi
          if(cmask(i) == 0 ) cycle

          if (read_daily_directly) then
            cloud(i) = xincldd(i) * 0.01
          else
c Here we add the NCEP cloud anomaly to the monthly anomaly from CRU05
c before converting percentage of cover to fraction
c 
            cloud(i) = (xincld(i,imonth) + xincldd(i)) * 0.01
          end if

c Bound cloud cover fraction between 0 and 1
          cloud(i) = min (cloud(i), 1.)
          cloud(i) = max (0.0, cloud(i))

        end do

        call derived_vars_from_cloud (jday)
        
      end if  ! read_radiation

      return
      end


c ---------------------------------------------------------------------
      subroutine create_daily_humidity (imonth)
c ---------------------------------------------------------------------
c
c This subroutine sets qd(i) (for i = 1 to npoi) when using real daily
c data
c
      use comgrid
      use compar
      use cominput
      use combcs
      use comatm
c
      implicit none
c
      include 'comsat.h'
c
c Arguments
c
      integer imonth
c
c Local variables
c
      integer i
      real 
     >  humidfrac,
     >  sphumid,
     >  qcov(12)

c Correction factors for covariances between relative humidity and
c temperature, so that one has the "true" monthly mean specific humidity
c (i.e., as would be calculated from the hourly specific humidity,
c rather than based on the monthly mean temperature). Correction factors
c for each month are based on data from the Trout Lake region of
c northern Wisconsin (for 1996-2000), where qcov = -1 +
c <RH*qs(Ta)>/<RH>*qs<Ta>.
      data qcov /0.194, 0.142, 0.102, -0.026, -0.017, 0.005, -0.006,
     >  -0.006, 0.017, 0.046, 0.052, 0.118/
 
c
      do i = 1, npoi
        if(cmask(i) == 0 ) cycle

        if (read_daily_directly) then
c Convert daily mean relative humidity to a fraction, then convert the
c fraction to daily mean specific humidity
          humidfrac = xinqd(i) / 100.
          qd(i) = qsat(humidfrac*esat(td(i)), psurf(i))
        else
c First convert relative humidity to a fraction and then convert the
c fraction to monthly mean specific humidity (based on the monthly mean
c temperature)
          humidfrac = xinq(i,imonth) / 100.
          sphumid = qsat(humidfrac*esat(273.16+xint(i,imonth)), psurf(i))

c Convert monthly mean specific humidity to daily mean using NCEP
c daily anomalies
          qd(i) = sphumid * xinqd(i)
        end if

        if (correct_sphumid) then
c Correct for covariances between relative humidity and temperature (see
c detailed comments above, regarding the qcov variable).
          qd(i) = qd(i) * (1. + qcov(imonth))
        end if

      end do

      return
      end
      

c ---------------------------------------------------------------------
      subroutine create_daily_wind (imonth)
c ---------------------------------------------------------------------
c
c This subroutine sets ud(i) (for i = 1 to npoi) when using real daily
c data
c
      use comgrid
      use compar
      use cominput
      use combcs
      use comatm
      use comveg
c
      implicit none
c
c Arguments
c
      integer imonth
c
c Local variables
c
      integer i
c 
      do i = 1, npoi
        if(cmask(i) == 0 ) cycle

        if (read_daily_directly) then
          ud(i) = xinwindd(i)
        else
c Here we multiply the NCEP fraction of windspeed by the CRU05
c climatological monthly windspeed
c 
          ud(i) = xinwind(i,imonth) * xinwindd(i)
        end if
c 

c WJS (02.10.10): Commented out the following adjustment because:
c  (1) It gave infinite wind speeds for crops (which have ztop(i,2)==0)
c  (2) There was supposed to be a later adjustment back down to a lower
c      height, but it looks like this readjustment never happened.
c So I got rid of this adjustment and instead changed siga back to
c 0.999, which is the approximate height of the wind speed data.
c
c$$$c TET change 10/29/02
c$$$c Input windspeed is assumed to be measured at canopy height
c$$$c (ztop(i,2))
c$$$c IBIS now has siga = 0.991 which is ~93m height.
c$$$c Adjust windspeed according to log-wind profile.
c$$$c So let
c$$$c displacement height, d = 0.65(ztop(i,2))
c$$$c roughness length, zm = 0.1(ztop(i,2))
c$$$c Equation found in Bonan (2002) and other texts:
c$$$c u(z2)-u(z1)=(frictionvelocity/k)*log((z2-d)/(z1-d))
c$$$c Use log-wind to solve for frictionvelocity and substitute
c$$$c into above equation:
c$$$c u(z2) = u(z1)*(1.+(log((z2-d)/(z1-d))) / (log((z1-d)/zm)))
c$$$c Use canopyheight = z1 = ztop(i,2), and z2=93m
c$$$c and substitute d and zm to get equation dependent on canopy height:
c$$$c
c$$$c WJS (01.21.10): Note that 1/(log((z1-d)/zm)) = 0.79824, regardless of
c$$$c the value of ztop(i,2); also, (z1 - d) = 0.35*ztop(i,2)
c$$$
c$$$        ud(i) = ud(i)*(1. + 0.79824*log((93.-0.65*ztop(i,2))/
c$$$     *    (0.35*ztop(i,2))))
c$$$c 
c$$$c Decided to use a canopy height of 20m over US after finding
c$$$c references of mature forest height between 10m and 35m (we are
c$$$c trying to simulate forests before human interference)
c$$$c with ztop(i,2) = 20m, the input windspeed is adjusted by a
c$$$c factor of 3.21
c$$$c
c$$$c WJS (01.21.10): As far as I can tell, this assumption of ztop = 20m
c$$$c doesn't come into play in the above equation
c 
      end do
      return
      end


c Other helper subroutines called by subroutine 'daily', which can be
c called either when using the weather generator or when using real
c daily data.

c ---------------------------------------------------------------------
      subroutine derived_vars_from_rads (jday)
c ---------------------------------------------------------------------
c
c If we already have values of rads(i) for every point i in the current
c day, then calculate today's values of trans(i) and cloud(i) for every
c point.
c
      use comgrid
      use compar
      use comatm
      use comwork
      use combcs,only: cmask
c
      implicit none
c
c Arguments
c
      integer jday     ! day of the year
c
c Local variables
c
      integer 
     >  i,             ! npoi loop index
     >  jj             ! latitude index of the current point

      real
     >  orbit,         ! earth's orbital angle (around the sun) in radians
     >  sw             ! effective solar constant
c
c Externals
c
      real calc_orbit, calc_sw  ! from toa_radiation.f


      orbit = calc_orbit(jday)
      sw = calc_sw(orbit)

      do i = 1, npoi
        if(cmask(i) == 0 ) cycle

c Calculate the solar transmission through the atmosphere using the
c ratio of measured daily solar radiation to theoretical
c top-of-atmosphere daily integrated solar radiation.
        jj = latindex(i)
        trans(i) = rads(i) / max(sw*coszen_mean(jj), epsilon)

c WJS 02.01.10: Constrain trans to be between 0.251 and 0.76, as it is
c when calculating trans from cloud cover (see derived_vars_from_clouds
c subroutine). Note that this constraint can mean that there is some
c inconsistency between the values of trans(i), rads(i), sw and
c coszen_mean(jj). One way to think about this is that there is some
c error in the estimation of sw (note that, as of 02.01.10, I don't
c think we use sw elsewhere in the code when read_radiation is true, so
c it's okay to mentally assign this error to sw). Alternatively, we
c could think of this constraint as applying not to trans per se but to
c the variables that are calculated from trans: namely, cloud (which is
c already constrained to be between 0 and 1) and fdiffuse.
        trans(i) = max(trans(i), 0.251)
        trans(i) = min(trans(i), 0.76)

c Estimate cloud cover based on trans, by simply using the inverse
c of the equation that we use below to estimate trans from cloud cover
c (in derived_vars_from_cloud).
        cloud(i) = min(1.0, max(0.0, 1. - ((trans(i) - 0.251) / 0.509)))
      end do

      return
      end
      

c ---------------------------------------------------------------------
      subroutine derived_vars_from_cloud (jday)
c ---------------------------------------------------------------------
c
c If we already have values of cloud(i) for every point i in the current
c day, then calculate today's values of trans(i) and rads(i) for every
c point.
c
c The code in this subroutine was moved from subroutine 'diurnal' 02.01.10
c
      use comgrid
      use compar
      use comatm
      use comwork
      use combcs,only: cmask
c
      implicit none
c
c Arguments
c
      integer jday     ! day of the year
c
c Local variables
c
      integer 
     >  i,             ! npoi loop index
     >  jj             ! latitude index of the current point

      real
     >  orbit,         ! earth's orbital angle (around the sun) in radians
     >  sw             ! effective solar constant
c
c Externals
c
      real calc_orbit, calc_sw  ! from toa_radiation.f


      orbit = calc_orbit(jday)
      sw = calc_sw(orbit)

      do i = 1, npoi
        if(cmask(i) == 0 ) cycle

c calculate the solar transmission through the atmosphere
c using simple linear function of tranmission and cloud cover
c 
c note that the 'cloud cover' data is typically obtained from
c sunshine hours -- not direct cloud observations
c 
c where, cloud cover = 1 - sunshine fraction 
c 
c different authors present different values for the slope and 
c intercept terms of this equation
c 
c Friend, A: Parameterization of a global daily weather generator for
c terrestrial ecosystem and biogeochemical modelling, Ecological 
c Modelling
c 
c Spitters et al., 1986: Separating the diffuse and direct component
c of global radiation and its implications for modeling canopy
c photosynthesis, Part I: Components of incoming radiation,
c Agricultural and Forest Meteorology, 38, 217-229.
c 
c A. Friend       : trans(i) = 0.251 + 0.509 * (1.0 - cloud(i))
c Spitters et al. : trans(i) = 0.200 + 0.560 * (1.0 - cloud(i))
c 
c we are using the values from A. Friend
c 
          trans(i) = 0.251 + 0.509 * (1.0 - cloud(i))
c
c calculate daily solar radiation at the surface based on daily
c integrated top-of-atmosphere solar radiation and transmission fraction
c
          jj = latindex(i)
          rads(i) = sw * coszen_mean(jj) * trans(i)
        end do

        return
        end


c ---------------------------------------------------------------------
c END HELPER SUBROUTINES CALLED BY SUBROUTINE 'DAILY'
c ---------------------------------------------------------------------



c
c ---------------------------------------------------------------------
      subroutine diurnal (istep, jday, plens, startp, endp, seed,
     >                    ilens, starti, endi)
c ---------------------------------------------------------------------
c
c uses:
c
      use comgrid
      use compar
      use cominput
      use comatm
      use comwork
      use comveg
      use comcrop
      use comnitr
      use combcs,only: cmask
c
      implicit none
c
c Arguments
c
      integer istep,   ! current time step (between 1 and niter)
     >        jday,    ! current day
     >        seed     !
c
      real
     >     plens,      ! length of the precip event (s)
     >     startp,     ! time to start precip event (s)
     >     endp,       ! time to end precip event (s)
     >     ilens,
     >     starti,
     >     endi
c
c determine the length of a precipitation event (between 4 and 24 hours),
c and time to start and end precipitation event. plen is in timesteps, while
c plens, startp, and endp are in seconds
c
c
c local variables
c
      integer i,         ! loop indice
     >        jj,
     >        ib         ! waveband number 1= visible, 2= near-IR
c
      real time,         ! time in seconds since midnight
     >     rtime,        ! time in hours
     >     orbit,        ! earth's orbital angle (around the sun) in radians
     >     xdecl,        ! solar declination angle
     >     xlat,         ! latitude in radians
     >     rads_now,     ! solar radiation at the surface, in this timestep (W m-2)
     >     fdiffuse,     ! fraction of indirect (diffuse) solar radiation
     >     fracw,         ! fraction of energy in each waveband
     >     gamma,        !
     >     ed,           ! daily average vapor pressure (n/m**2)
     >     tdewd,        ! daily average dew point temperature (K)
     >     tdew_min,     ! minimum daily dew point temperature (K)
     >     tdew_max,     ! maximum daily dew point temperature (K)
     >     tdew_now,     ! current dew point temperature (K)
     >     qmin,         !
     >     qmax,
     >     qsa,
     >     e_cur,           ! vapor pressure in this time step (n/m**2)
     >     ran2,
     >     emb,
     >     ea,
     >     ec,
     >     dtair,
     >     dtcloud,
     >     mycloud       ! cloud cover for the current point, either measured or estimated
c
      integer checkP,
     >        niter,
     >        plen,
     >        plenmin,
     >        plenmax
c
      include 'comsat.h'
c
c Externals
c
      real get_time  ! from utilities.f
      real calc_orbit, calc_xdecl  ! from toa_radiation.f
      real tdew_from_e, e_from_tdew  ! from humid.f

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
c
 9001 continue
c
c do for all gridcells
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c ---------------------------------------------------------------------- 
c * * * solar calculations * * *
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
     >    (acos(-(sin(xlat)*sin(xdecl)) / (cos(xlat)*cos(xdecl))))
c
c WJS 02.01.10: Calculation of 'trans' moved to
c 'derived_vars_from_cloud', called by subroutine 'daily'.
c
c Estimate solar radiation at the surface in the current timestep using
c daily radiation and the ratio of (top-of-atmosphere radiation in this
c timestep) to daily-integrated top-of-atmosphere radiation.
c
c top-of-atmosphere radiation in this timestep = sw*coszen(i)
c daily-integrated top-of-atmosphere radiation = sw*coszen_mean(jj)
c
        rads_now = rads(i) * coszen(i) / max(coszen_mean(jj), epsilon)
c
c calculate the fraction of indirect (diffuse) solar radiation
c based upon the cloud cover
c
c WJS 01.28.10: Now that we sometimes calculate trans based on measured
c solar radiation, rather than based on cloud cover, I think the above
c comment should say, "calculate the fraction of indirect (diffuse)
c solar radiation based upon the solar transmission fraction (trans)"
c
c note that these relationships typically are measured for either
c monthly or daily timescales, and may not be exactly appropriate
c for hourly calculations -- however, in ibis, cloud cover is fixed
c through the entire day so it may not make much difference
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
        if (trans(i).gt.0.75) fdiffuse = 0.166
    
c in case trans is near 0, constrain fdiffuse to be <= 1
c note that, with the above equation, fdiffuse can never be < 0, so we
c don't need to check for that.
        fdiffuse = min (fdiffuse, 1.)
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
        do 120 ib = 1, nband
c
c calculate the fraction in each waveband
c
          fracw = 0.46 + 0.08 * float(ib - 1)
c
c calculate the direct and indirect solar radiation
c
          solad(i,ib) = fracw * rads_now * (1. - fdiffuse)
c
          solai(i,ib) = fracw * rads_now * fdiffuse
c
  120   continue
c
c ---------------------------------------------------------------------- 
c * * * temperature calculations * * *
c ---------------------------------------------------------------------- 
c
c assign hourly temperatures using tmax and tmin 
c following Environmental Biophysics, by Campbell and Norman, p.23
c
c WJS (4-27-10): But note that the value of tmax_wt that we now use is different from the value in C&N
c
c this function fits a fourier series to the diurnal temperature cycle
c note that the maximum temperature occurs at 2:00 pm local solar time
c
c note that the daily mean value of gamma is given by tmax_wt
c
        gamma = tmax_wt - 0.46 * sin (      pi / 12.0 * rtime + 0.9) 
     >               + 0.11 * sin (2.0 * pi / 12.0 * rtime + 0.9)
c
        ta(i) = tmax(i) * gamma + tmin(i) * (1.0 - gamma)
c
c ---------------------------------------------------------------------- 
c * * * humidity calculations * * *
c ---------------------------------------------------------------------- 
c
        if (conserve_tdew) then  ! conserve daily mean dew point temperature
c 
c adjust dew point temperatures against daily minimum temperatures
c
c To do this, dew point temperature is written as an approximate sine
c function (same as ta) to preserve the daily mean dew point
c temperature, while also preventing RH from exceeding 99% at night. If
c daily average dew point temperature is lower than tmin (actually,
c slightly lower because of the assumption that RH cannot exceed 99%),
c then we assume constant humidity throughout the day; but if daily
c average dew point temperature is greater than tmin, then we impose a
c sinusoidal cycle on humidity that ensures that RH doesn't exceed 99%
c at night.
c
c Note that the daily mean RH is *not* preserved, and therefore the
c output RH will be slightly different from what was read in
c 
c first calculate daily average vapor pressure (ed) from qd(i), then
c convert this to daily average dew point temperature (tdewd)
c
          ed = e_from_q (qd(i), psurf(i))
          tdewd = tdew_from_e (ed)
c
c calculate the dew point temperature corresponding to 99% RH when the
c temperature is tmin(i)
c
          tdew_min = tdew_from_e (0.99 * esat(tmin(i)))
c
c But, if daily average dew point temperature is already lower than the
c above-calculated tdew_min, then use daily average dew point
c temperature as the min. daily dew point temperature; in this case,
c dew point temperature (and thus specific humidity) will be constant
c throughout the day. 
c
          tdew_min = min(tdewd, tdew_min)
c
c Calculate maximum daily tdew assuming the same sinusoidal curve as is
c used  for temperature. Note that if tdew_min = tdewd, then tdew_max
c will also equal tdewd.
c
          tdew_max = (tdewd - tmin_wt * tdew_min) / tmax_wt
c
c Calculate the hourly dew point temperature
c
          tdew_now = tdew_max * gamma + tdew_min * (1.0 - gamma)
c
c Finally, calculate the hourly specific humidity corresponding to this
c dew point temperature (note: although we are using a function called
c 'qsat', this function works equally well on non-saturation values of
c specific humidity, as is done here).
c
          qa(i) = qsat (e_from_tdew (tdew_now), psurf(i))

        else  ! .not. conserve_tdew: conserve daily mean specific humidity

c adjust specific humidity against daily minimum temperatures
c
c To do this, qa is written as an approximate sine function (same as ta)
c to preserve the daily mean specific humidity, while also preventing rh
c from exceeding 99% at night
c
c Note that the daily mean RH is *not* preserved, and therefore the
c output RH will be slightly different from what was read in.
c
c first adjust so that maximum RH cannot exceed 99% at night
c
          qmin = min (qd(i), qsat(0.99*esat(tmin(i)), psurf(i)))
          qmax = (qd(i) - tmin_wt * qmin) / tmax_wt
c
c calculate the hourly specific humidity, using the above adjustments
c
          qa(i) = qmax * gamma + qmin * (1.0 - gamma)

        end if
c 
c if needed, adjust humidity to no greater than 99% at other times of
c the day (in which case the daily mean *specific* humidity is also not
c preserved)
c
        qsa  = qsat(0.99*esat(ta(i)), psurf(i))
        qa(i) = min (qsa, qa(i))
c
c calculate the hourly relative humidity 
c
        e_cur = e_from_q (qa(i), psurf(i))
        rh(i) = 100.0 * e_cur / esat(ta(i))
c
c ---------------------------------------------------------------------- 
c * * * wind speed calculations * * *
c ---------------------------------------------------------------------- 
c
c following logic of the EPIC weather generator
c select random wind speed following this equation
c
        ua(i) = 1.13989 * ud(i) * (-log(ran2(seed)))**0.30 
c
c fix wind speeds to always be above 2.5 m/sec and below 25.0 m/sec
c
c WJS 02.10.10: The original constraint was 2.5 - 10; based on
c discussions with CJK, we think 10 is too low, so we're using a max of
c 25 (somewhat arbitrarily)
c
c        ua(i) = max (2.5, min (10.0, ua(i)))
        ua(i) = max (2.5, min (25.0, ua(i)))
c
c ---------------------------------------------------------------------- 
c * * * ir flux calculations * * *
c ---------------------------------------------------------------------- 
c
c clear-sky emissivity as a function of water vapor pressure
c and atmospheric temperature
c
c calculate the ir emissivity of the clear sky
c using equation from idso (1981) water resources res., 17, 295-304
c
        emb = 0.01 * (psurf(i) * qa(i) / (0.622 + qa(i)))
        ea  = 0.70 + 5.95e-5 * emb * exp (1500.0 / ta(i))
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
        fira(i) = (1. -  cloud(i)) * ea * stef * (ta(i) - dtair  )**4 +
     >                   cloud(i)  * ec * stef * (ta(i) - dtcloud)**4
c
c ---------------------------------------------------------------------- 
c * * * snow and rain calculations * * *
c ---------------------------------------------------------------------- 
c
c reset snow and rain to zero
c
        snowa(i) = 0.0
        raina(i) = 0.0
c
c determine the number of timesteps per day
c
        niter = int (86400.0 / dtime)
c
c change the rain length when the amount of rainfall/timestep is
C too high (at the first time step)
c
        if (time.lt.dtime) then
c
           plen = plens / dtime
           plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
           plenmax = max (int (24.0 * 3600. / dtime), plenmin)
           checkP = 0
c
           do  while (((precip(i)/plen) .gt. 15).and.(plen.lt.plenmax))
              plen = plen + 1
              checkP = 1
           end do
c
           !if (checkP.eq.1) then
           if ((checkP.eq.1) .and. .not.(error_showed_plen)) then   !Y.Li
c
              print *, 'WARNING: plen changed', i,
     $             int(precip(i)), int(plens/dtime), plen
              plens = dtime * plen
              startp = dtime * min (niter-plen,
     >             int(ran2(seed)*(niter-plen+1)))
              endp = startp + plen *dtime

              error_showed_plen = .true.
              goto 9001
           end if
c
        end if
c
c if precipitation event then calculate
c
        if (time.ge.startp .and. time.lt.endp) then  
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
          if (ta(i)-273.15 .gt. 2.5) then
            raina(i) = precip(i) / plens
          else
            snowa(i) = precip(i) / plens
          endif
c
        endif

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
c
        if (time.ge.starti .and. time.lt.endi
     >      .and. irrigate(i) .eq. 1
     >      .and. precip(i) .eq. 0.00) then  
c
          xirriga(i) = xirrig(i) / ilens
c
c update annual total - totirrig
c rate of irrigation multiplied by length of timestep (mm/s * s) = total applied
c for this timestep 
c
          totirrig(i) = totirrig(i) + (xirriga(i) * dtime)
c
        endif 
c
  100 continue
c
      return
      end
c

c -----------------------------------------------------------------------------------
      subroutine dailymet (imonth, iday, seed, jdaily)
c -----------------------------------------------------------------------------------
c
c overview
c
c this routine generates daily weather conditions from monthly-mean
c climatic parameters
c
c specifically, this routine generates daily values of
c
c  - daily total precipitation
c  - daily maximum temperature
c  - daily minimum temperature
c  - daily average cloud cover
c  - daily average relative humidity
c  - daily average wind speed
c
c in order to generate daily weather conditions, the model uses a series
c of 'weather generator' approaches, which generate random combinations of
c weather conditions based upon the climatological conditions
c
c in general, this weather generator is based upon the so-called Richardson
c weather generator
c
c appropriate references include:
c
c Geng, S., F.W.T. Penning de Vries, and L. Supit, 1985:  A simple
c method for generating rainfall data, Agricultural and Forest
c Meteorology, 36, 363-376.
c
c Richardson, C. W. and Wright, D. A., 1984: WGEN: A model for 
c generating daily weather variables: U. S. Department of
c Agriculture, Agricultural Research Service.
c
c Richardson, C., 1981: Stochastic simulation of daily
c precipitation, temperature, and solar radiation. Water Resources 
c Research 17, 182-190.
c
c uses:
c
      use comgrid
      use compar
      use combcs
      use comatm
      use comsum
      use comveg
      use comsoi
      use comcrop
c
      implicit none
c
c Arguments
c
      integer seed, 
     >        iyear,
     >        imonth,
     >        iday,
     >        jdaily   ! 1 if reading in daily weather data
     >                 ! 0 if using random/statistical weather generator      
c
c local variables
c
      integer it1w,      ! indice of previous month (interpolation)
     >        it2w,      ! indice of following month (interpolation)
     >        i,j, k     ! loop indice
c
      real rwork,          ! 
     >     omcloud,        ! cloud cover
     >     omqd,           ! humidity
     >     omtmax,         ! maximum temperature
     >     ran2,           ! function random number generator
     >     dt,             ! used for interpolation
     >     pwet,           ! monthly-average probability of rainy day
     >     pwd,            ! probability of a wet day after a dry day
     >     pww,            ! probability of a wet day after a wet day
     >     rndnum,         ! random number to decide if wet or dry day
     >     rainpwd,        ! average rainfall per wet day
     >     alpha,          ! parameter for gamma function
     >     beta,           ! parameter for gamma function
     >     aa,
     >     ab,
     >     tr1,
     >     tr2,
     >     rn1,rn2, rn3,rn,  !random numbers
     >     s1,
     >     s2,
     >     s12,
     >     z,
     >     tdm,            ! mean daily mean temperature
     >     trngm,          ! mean daily mean temperature
     >     tmaxm,          ! mean maximum temperature
     >     tminm,          ! mean minimum temperature
     >     tmaxd,          ! maximum temperatures for dry days
     >     tmaxw,          ! maximum temperatures for wet days
     >     tmaxe,          !'expected' maximum temperature for today
     >     tmine,          !'expected' minimum temperature for today
     >     tmaxs,          ! standard deviation in maximum temperature (K)
     >     tmins,          ! standard deviation in minimum temperature (K)
     >     cloudm,         ! mean cloud cover for today (fraction)
     >     cloudd,         ! dry day cloud cover
     >     cloudw,         ! wet day cloud cover
     >     cloude,         ! expected cloud cover today
     >     clouds,         ! standard deviation of cloud fraction
     >     v,
     >     tdum,           ! storage variable
     >     qdm,            ! mean relative humidity
     >     qdd,            ! dry day relative humidity
     >     qdw,            ! wet day relative humidity 
     >     qde,            ! expected relative humidity (based on wet/dry decision)
     >     qdup,           ! upper bound of humidity distribution function
     >     qdlow,          ! lower bound of humidity distribution function
     >     y,
     >     b3,
     >     b2,
     >     b1,
     >     x1,
     >     amn,
     >     eud             ! expected daily average wind speed from monthly mean
 
c
      real a(3,3),
     >     b(3,3)
c
      real
     >     ee(3), 
     >     r(3),
     >     rr(3),
     >     x(3)
c
      integer
     >     iyrlast,
     >     nrun
c
      real
     >     precipfac,
     >     dif,
     >     humidfrac,
     >     sphumid
c
c Externals
c
      real gdd_natveg, gdd_crop, gdd_soil  ! from gdd.f
c
c define autocorrelation matrices for Richardson generator
c
c note that this matrix should be based upon a statistical
c analysis of regional weather patterns
c
c for global simulations, we use 'nominal' values
c
      data a /  0.600,  0.500,  0.005,
     >          0.010,  0.250,  0.005, 
     >          0.020,  0.125,  0.250 /
c
      data b /  0.500,  0.250, -0.250,
     >          0.000,  0.500,  0.250, 
     >          0.000,  0.000,  0.500 /
c
      include 'comsat.h'
c
c ---------------------------------------------------------------------- 
c * * * initial setup for daily climate calculations * * *
c ---------------------------------------------------------------------- 
c
c define working variables
c
      rwork = (grav / rair / 0.0065)
c
c 'omega' parameters used to calculate differences in expected
c climatic parameters on wet and dry days
c
c following logic of weather generator used in the EPIC crop model
c
c omcloud -- cloud cover
c omqd    -- humidity
c omtmax  -- maximum temperature
c
      omcloud = 0.90    ! originally 0.90
      omqd    = 0.50    ! originally 0.50
      omtmax  = 0.75    ! originally 0.75
c
c calculate weighting factors used in interpolating climatological
c monthly-mean input values to daily-mean values
c
c this is a simple linear interpolation technique that takes into
c account the length of each month
c
c      if (jdaily .eq. 0) then
c        if (float(iday).lt.float(ndaypm(imonth)+1)/2.0) then
c          it1w = imonth - 1
c          it2w = imonth 
c          dt   = (float(iday) - 0.5) / ndaypm(imonth) + 0.5
c        else
c         it1w = imonth 
c          it2w = imonth + 1
c          dt   = (float(iday) - 0.5) / ndaypm(imonth) - 0.5
c        end if
c
c        if (it1w.lt. 1) it1w = 12
c        if (it2w.gt.12) it2w = 1
c      end if
c
c initialize this year's values of gdd0, gdd5, tc, tw
c
      if (iday.eq.1 .and. imonth.eq.1) then
c
        call const (tcthis, npoi,  100.0)
        call const (twthis, npoi, -100.0)
c
        call const (gdd0this, npoi, 0.0)
        call const (gdd5this, npoi, 0.0)
        call const (gdd0cthis, npoi, 0.0)
        call const (gdd8this, npoi, 0.0)
        call const (gdd10this, npoi, 0.0)
c
c constant for gdd based on 10 cm soil temperature (planting bed-
c used to calculate leaf emergence after planting
c
c for all crop plant functional types
c and aboveground annual npp for crops
c 
c constant for gdd based on 10 cm soil temperature (planting bed-
c used to calculate leaf emergence after planting
c
        do 20 i = 1, npoi
          if(cmask(i) == 0 ) cycle

          do 50 j = 1, npft
c
            if (j .le. scpft-1) then    ! natural vegetation  
              ayanpp(i,j)     = 0.0
            else if (croppresent(i,j) .eq. 0 .and. j .ge. scpft) then
              gddplant(i,j)   = 0.0
              gddtsoi(i,j)    = 0.0
              ayanpp(i,j)     = 0.0
            endif
 50       continue
 20     continue
c           
      endif
c
c ---------------------------------------------------------------------- 
c * * * set daily climatic variables for entire domain * * *
c ---------------------------------------------------------------------- 
c
      do 200 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c 
c ---------------------------------------------------------------------- 
c * * * other daily climate calculations * * *
c ---------------------------------------------------------------------- 
c
c calculated temperature extremes -- for vegetation limits (deg c)
c
c for this purpose, use the 10-day running mean temperature
c
        tcthis(i) = min (tcthis(i), (a10td(i) - 273.16))
        twthis(i) = max (twthis(i), (a10td(i) - 273.16))
c
c update this year's growing degree days
c
        gdd0this(i)  = gdd0this(i)  + gdd_natveg(td(i)-c_to_k, 0.0)
        gdd0cthis(i) = gdd0cthis(i) + gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, 0.0, 26.0)  ! wheat
        gdd5this(i)  = gdd5this(i)  + gdd_natveg(td(i)-c_to_k, 5.0)
        gdd8this(i)  = gdd8this(i)  + gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, 8.0, 30.0)  ! maize (prior to 4-22-10)
        gdd10this(i) = gdd10this(i) + gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, 10.0, 30.0) ! soybean (& maize, since 4-22-10)
c
c accumulate growing degree days for planted crops past planting
c
        do 150 j = scpft, ecpft 
c
c for crops except winter wheat
c 
c We use croppresent rather than croplive, so we continue to accumulate GDD between maturity & harvest,
c in case any logic depends on that continued accumulation of GDD
c
          if (croppresent(i,j) .eq. 1 .and. iwheat .ne. 2) then
c
            gddplant(i,j) = gddplant(i,j) + 
     >        gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c
            gddtsoi(i,j)  = gddtsoi(i,j) + 
     >        gdd_soil(tsoi(i,1)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c
c if winter wheat is planted, reduce thermal time accumulation
c by vernalization factor (calculated in crops.f) until crop
c is fully vernalized (cold-hardened)
c
c We use croppresent rather than croplive, so we continue to accumulate GDD between maturity & harvest,
c in case any logic depends on that continued accumulation of GDD
c
          else if (croppresent(i,j) .eq. 1 .and. iwheat .eq. 2) then
c
            gddplant(i,j) = gddplant(i,j) + 
     >        vf(i) * gdd_crop(td(i)-c_to_k, tmax(i)-c_to_k, tmin(i)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c
            gddtsoi(i,j)  = gddtsoi(i,j) + 
     >        vf(i) * gdd_soil(tsoi(i,1)-c_to_k, baset(j)-c_to_k, mxtmp(j))
c
        endif
c
 150  continue
c
 200  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine diurnalmet (istep, jday, plens, startp, endp, seed,
     >                       ilens, starti, endi)
c ---------------------------------------------------------------------
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comveg
      use comwork
      use comhour
      use comcrop
      use comnitr
      use combcs,only: cmask
c
      implicit none
c
c Arguments
c
      integer istep,   ! current time step (between 1 and niter)
     >        jday,    ! current day
     >        seed     !
c
      real
     >     plens,      ! length of the precip event (s)
     >     startp,     ! time to start precip event (s)
     >     endp,       ! time to end precip event (s)
     >     ilens,
     >     starti,
     >     endi,
     >     truecloud 
c
      integer i,         ! loop indice
     >        jj, 
     >        ib         ! waveband number 1= visible, 2= near-IR
c
      real time,         ! time in seconds since midnight
     >     rtime,        ! time in hours
     >     orbit,        ! earth's orbital angle (around the sun) in radians
     >     xdecl,        ! solar declination angle
     >     sw,           ! effective solar constant
     >     xlat,         ! latitude in radians
     >     fdiffuse,     ! fraction of indirect (diffuse) solar radiation
     >     fracw,        ! fraction of energy in each waveband
     >     wfrac,
     >     gamma,        !
     >     qmin,         !
     >     qmax,
     >     qsa,
     >     e_cur,           ! vapor pressure in this time step (n/m**2)
     >     ran2,
     >     emb,
     >     ea,
     >     ec,
     >     dtair,
     >     dtcloud
c
      integer checkP,
     >        niter,
     >        plen,
     >        plenmin,
     >        plenmax
c
      include 'comsat.h'
c
c Externals
c
      real get_time  ! from utilities.f
      real calc_orbit, calc_xdecl, calc_sw  ! from toa_radiation.f
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
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c ---------------------------------------------------------------------- 
c * * * solar calculations * * *
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
     >    (acos(-(sin(xlat)*sin(xdecl)) / (cos(xlat)*cos(xdecl))))
c
c calculate the solar transmission through the atmosphere
c using simple linear function of tranmission and cloud cover
c
c note that the 'cloud cover' data is typically obtained from
c sunshine hours -- not direct cloud observations
c
c where, cloud cover = 1 - sunshine fraction 
c
c different authors present different values for the slope and 
c intercept terms of this equation
c
c Friend, A: Parameterization of a global daily weather generator for
c terrestrial ecosystem and biogeochemical modelling, Ecological 
c Modelling
c
c Spitters et al., 1986: Separating the diffuse and direct component
c of global radiation and its implications for modeling canopy
c photosynthesis, Part I: Components of incoming radiation,
c Agricultural and Forest Meteorology, 38, 217-229.
c
c A. Friend       : trans(i) = 0.251 + 0.509 * (1.0 - cloud(i))
c Spitters et al. : trans(i) = 0.200 + 0.560 * (1.0 - cloud(i))
c
c we are using the values from A. Friend
c
cjk     trans(i) = 0.251 + 0.509 * (1.0 - cloud(i)) 
        trans(i) = cloud(i) / sw     ! cloud(i) is surface insolation
c
cjk
c        if (jday .eq. 1) open (26, file='met.trans',status='unknown')
c        write(26,*) jday, istep, trans(i)
c
c calculate the fraction of indirect (diffuse) solar radiation
c based upon the cloud cover
c
c note that these relationships typically are measured for either
c monthly or daily timescales, and may not be exactly appropriate
c for hourly calculations -- however, in ibis, cloud cover is fixed
c through the entire day so it may not make much difference
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
        if (trans(i).gt.0.75) fdiffuse = 0.166
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
        do 120 ib = 1, nband
c
c calculate the fraction in each waveband
c
          wfrac = 0.46 + 0.08 * float(ib - 1)
c
c calculate the direct and indirect solar radiation
c
cjk       solad(i,ib) = sw * coszen(i) * wfrac * trans(i) *
cjk  >                  (1. - fdiffuse)  
cjk
          solad(i,ib) = wfrac * cloud(i) * (1. - fdiffuse) 
c
          solai(i,ib) = wfrac * cloud(i) * fdiffuse
c
  120   continue
c
c ---------------------------------------------------------------------- 
c * * * temperature calculations * * *
c ---------------------------------------------------------------------- 
c
c assign hourly temperatures using tmax and tmin 
c following Environmental Biophysics, by Campbell and Norman, p.23
c
c WJS (4-27-10): But note that the value of tmax_wt that we now use is different from the value in C&N
c
c this function fits a fourier series to the diurnal temperature cycle
c note that the maximum temperature occurs at 2:00 pm local solar time
c
c note that the daily mean value of gamma is given by tmax_wt
c
        gamma = tmax_wt - 0.46 * sin (      pi / 12.0 * rtime + 0.9) 
     >               + 0.11 * sin (2.0 * pi / 12.0 * rtime + 0.9)
c
cjk    ta(i) = tmax(i) * gamma + tmin(i) * (1.0 - gamma)
c
c ---------------------------------------------------------------------- 
c * * * humidity calculations * * *
c ---------------------------------------------------------------------- 
c
c adjust specific humidity against daily minimum temperatures
c
c To do this, qa is written as an approximate sine function (same as ta)
c to preserve the daily mean specific humidity, while also preventing rh
c from exceeding 99% at night
c
c Note that the daily mean RH is *not* preserved, and therefore the
c output RH will be slightly different from what was read in.
c
c first adjust so that maximum RH cannot exceed 99% at night
c
cjk     qmin = min (qd(i), qsat(0.99*esat(tmin(i)), psurf(i)))
cjk     qmax = (qd(i) - tmin_wt * qmin) / tmax_wt
c
c if needed, adjust again to 99% at other times of the day (in which
c case the daily mean *specific* humidity is also not preserved)
c
        qsa  = qsat(0.99*esat(ta(i)), psurf(i))
c
c calculate the hourly specific humidity, using the above adjustments
c
cjk     qa(i) = min (qsa, qmax * gamma + qmin * (1.0 - gamma))
        qa(i) = min (qsa, qd(i))
c
c calculate the hourly relative humidity 
c
        e_cur = e_from_q (qa(i), psurf(i))
        rh(i) = 100.0 * e_cur / esat(ta(i))
c
c ---------------------------------------------------------------------- 
c * * * wind speed calculations * * *
c ---------------------------------------------------------------------- 
c
c following logic of the EPIC weather generator
c select random wind speed following this equation
c
cjk     ua(i) = 1.13989 * ud(i) * (-log(ran2(seed)))**0.30 
c
c fix wind speeds to always be above 2.5 m/sec and below 25.0 m/sec
c
c WJS 02.10.10: The original constraint was 2.5 - 10; based on
c discussions with CJK, we think 10 is too low, so we're using a max of
c 25 (somewhat arbitrarily)
c
c        ua(i) = max (2.5, min (10.0, ua(i)))
        ua(i) = max (2.5, min (25.0, ua(i)))
c
c ---------------------------------------------------------------------- 
c * * * ir flux calculations * * *
c ---------------------------------------------------------------------- 
c
c clear-sky emissivity as a function of water vapor pressure
c and atmospheric temperature
c
c calculate the ir emissivity of the clear sky
c using equation from idso (1981) water resources res., 17, 295-304
c
        emb = 0.01 * (psurf(i) * qa(i) / (0.622 + qa(i)))
        ea  = 0.70 + 5.95e-5 * emb * exp (1500.0 / ta(i))
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
c cjk 
c
        truecloud = 1. - ((trans(i) - 0.251) / 0.509) 
        fira(i) = (1. -  truecloud) * ea * stef * (ta(i) - dtair  )**4 +
     >                   truecloud  * ec * stef * (ta(i) - dtcloud)**4
c
cc      fira(i) = (1. -  cloud(i)) * ea * stef * (ta(i) - dtair  )**4 +
cc   >                   cloud(i)  * ec * stef * (ta(i) - dtcloud)**4
c
c ---------------------------------------------------------------------- 
c * * * snow and rain calculations * * *
c ---------------------------------------------------------------------- 
c
c reset snow and rain to zero
c
        snowa(i) = 0.0
        raina(i) = 0.0
c
c determine the number of timesteps per day
c
        niter = int (86400.0 / dtime)
c
c change the rain length when the amount of rainfall/timestep is
C too high (at the first time step)
c
c        if (time.lt.dtime) then
c
c           plen = plens / dtime
c           plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
c           plenmax = max (int (24.0 * 3600. / dtime), plenmin)
c           checkP = 0
c
c           do  while (((precip(i)/plen) .gt. 15).and.(plen.lt.plenmax))
c              plen = plen + 1
c              checkP = 1
c           end do
c
c           if (checkP.eq.1) then
c
c              print *, 'WARNING: plen changed', i,
c     $             int(precip(i)), int(plens/dtime), plen
c              plens = dtime * plen
c              startp = dtime * min (niter-plen,
c     >             int(ran2(seed)*(niter-plen+1)))
c              endp = startp + plen *dtime
c              goto 9001
c           end if
c
c        end if
c
c if precipitation event then calculate
c
c cjk        if (time.ge.startp .and. time.lt.endp) then  
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
          if (ta(i)-273.15 .gt. 2.5) then
            raina(i) = precip(i) / plens 
          else
            snowa(i) = precip(i) / plens
          endif
c
c cjk        endif
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
c
        if (time.ge.starti .and. time.lt.endi
     >      .and. irrigate(i) .eq. 1
     >      .and. precip(i) .eq. 0.00) then  
c
          xirriga(i) = xirrig(i) / ilens
c
c update annual total - totirrig
c rate of irrigation multiplied by length of timestep (mm/s * s) = total applied
c for this timestep 
c
          totirrig(i) = totirrig(i) + (xirriga(i) * dtime)
c
        endif 
c
  100 continue
c
      return
      end
c
