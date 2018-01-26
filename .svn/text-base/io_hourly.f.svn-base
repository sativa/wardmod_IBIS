c io_hourly.f
c Created by Jason Patton 2010-06-22
c
c Holding place for subprograms dealing with outputting hourly data

c ----------------------------------------
      SUBROUTINE openHourlyOut
c ----------------------------------------
c Created by Jason Patton 2010-06-22
c
c Opens files for output
c Requires directory output/hourly to exist

      OPEN(1277, FILE = 'output/hourly/crop_properties.out', 
     >     STATUS = 'REPLACE')
      OPEN(1278, FILE = 'output/hourly/soil_properties.out', 
     >     STATUS = 'REPLACE')
c      OPEN(1279, FILE = 'flux_properties.out', STATUS = 'REPLACE')
c      OPEN(1280, FILE = 'canopy_properties.out', STATUS = 'REPLACE')
      OPEN(1281, file = 'output/hourly/gdd_properties.out',
     >     STATUS = 'REPLACE')

      END SUBROUTINE openHourlyOut
c ----------------------------------------

c ----------------------------------------
      SUBROUTINE closeHourlyOut
c ----------------------------------------
c Created by Jason Patton 2010-06-22
c
c Closes files for output

      CLOSE(1277)
      CLOSE(1278)
c      CLOSE(1279)
c      CLOSE(1280)
      CLOSE(1281)

      END SUBROUTINE closeHourlyOut
c ----------------------------------------


c ----------------------------------------
      SUBROUTINE whourly (thisYear, thisMonth, thisDay, thisStep)
c ----------------------------------------
c Created by Jason Patton 2010-06-22
c
c Writes output

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
      USE comsoi
      USE comveg
      USE com1d

      IMPLICIT NONE

      INTEGER thisYear, thisMonth, thisDay, thisStep
      INTEGER thisLayer
      INTEGER thisPFT
      INTEGER thisPoint

c Determine PFTs (only corn or soybean right now)
      IF(exist(1,13) .GE. 1) THEN
         thisPFT = 13
      ELSEIF(exist(1,14) .GE. 1) THEN
         thisPFT = 14
      ELSE
         PRINT*, 'no crops!'
      ENDIF

c Determine grid point to use (grabs middle point of domain)
      thisPoint = npoi/2

c Check for end of day for things that are only output on daily timestep

      IF(thisStep .EQ. 24) THEN

c Output:
c Biomass, daily timesteps, biomass(npoi, npft)
c Leaf area index, daily timesteps plai(npoi, npft)
c Canopy height, daily timesteps, ztop(npoi,canopy), canopy = 1 for crops

         WRITE(1277, '(I4, 2(1X, I2.2), 3(1X, F12.5))')
     >        thisYear, thisMonth, thisDay, ! crop_properties
     >        biomass(thisPoint, thisPFT),  ! biomass (kg m-2, carbon only)
     >        plai(thisPoint, thisPFT),     ! lai (m3 m-3)
     >        ztop(thisPoint, 1)            ! canopy height (m)

c Output:
c GDD (plant), daily timesteps, gddplant(npoi, npft)
c GDD (soil), daily timesteps, gddsoi(npoi, npft)

         WRITE(1281, '(I4, 2(1X, I2.2), 4(1X, F12.5))')
     >        thisYear, thisMonth, thisDay, ! gdd_properties
     >        gddplant(thisPoint, thisPFT), ! gdd (plant, C) 
     >        gddtsoi(thisPoint, thisPFT)   ! gdd (soil, C)

      ENDIF

c Output:
c Soil moisture, hourly timesteps, wsoi(npoi, nsoilay)
c Soil temperature, hourly timesteps, tsoi(npoi, nsoilay)

      WRITE(1278, '(I4, 3(1X, I2.2))', ADVANCE = 'NO')
     >     thisYear, thisMonth, thisDay, (thisStep - 1) ! soil_properties

c Loop over all soil layers and add to output
      DO thisLayer = 1, nsoilay
         WRITE(1278, '(2(1X, F12.5))', ADVANCE = 'NO')
     >        wsoi(thisPoint, thisLayer), ! soil moisture (pore space fraction)
     >        tsoi(thisPoint, thisLayer)  ! soil temperature (K)
      ENDDO

      WRITE(1278, *)

c     Albedo, hourly timesteps, asurd(npoi, nband), asuri(npoi, nband)
c     IR, SH, LH, etc. look at com1d.f
c     CO2 flux, hourly timesteps, tneetot(npoi), tco2mic(npoi), tco2root(npoi)

c     In-canopy temperature, hourly timesteps, t34(npoi)
c     In-canopy spec. humidity, hourly timesteps, q34(npoi)
c     Dew (actually, evaporation from leaves and stems), hourly timesteps, fvaplw(npoi)

      END SUBROUTINE whourly
