c This file contains functions & subroutines for calculating top-of
c -atmosphere radiation and related variables.
c
c Created by Bill Sacks, 01.29.10; this code was extracted from
c weather.f : diurnal

c ---------------------------------------------------------------------
      real function calc_orbit (jday)
c ---------------------------------------------------------------------
c
c Return the earth's orbital angle (around the sun) in radians
c
      use compar        ! just needed for pi
c
      implicit none
c
c Arguments
c
      integer jday    ! current day

      calc_orbit = 2.0 * pi * float(jday) / 365.2425
      
      return
      end

c ---------------------------------------------------------------------
      real function calc_angle (rtime)
c ---------------------------------------------------------------------
c
c Return the solar hour angle in radians
c
      use compar        ! just needed for pi
c
      implicit none
c
c Arguments
c
      real rtime   ! time in hours

      calc_angle = 2.0 * pi * (rtime - 12.0) / 24.0

      return
      end

c ---------------------------------------------------------------------
      real function calc_xdecl (orbit)
c ---------------------------------------------------------------------
c
c Calculate the current solar declination angle
c Ref: Global Physical Climatology, Hartmann, Appendix A
c
      implicit none
c
c Arguments
c 
      real orbit   ! earth's orbital angle (around the sun) in radians

      calc_xdecl =  0.006918                    -
     >              0.399912 * cos(orbit)       + 
     >              0.070257 * sin(orbit)       -
     >              0.006758 * cos(2.0 * orbit) +
     >              0.000907 * sin(2.0 * orbit) -
     >              0.002697 * cos(3.0 * orbit) +
     >              0.001480 * sin(3.0 * orbit)

      return
      end

c ---------------------------------------------------------------------
      real function calc_sw (orbit)
c ---------------------------------------------------------------------
c
c Calculate the effective solar constant, including effects of eccentricity
c Ref: Global Physical Climatology, Hartmann, Appendix A
c
      implicit none
c
c Arguments
c 
      real orbit   ! earth's orbital angle (around the sun) in radians

      calc_sw = 1370. * (1.000110 +
     >                   0.034221 * cos(orbit)        + 
     >                   0.001280 * sin(orbit)        + 
     >                   0.000719 * cos(2.0 * orbit)  +
     >                   0.000077 * sin(2.0 * orbit)) 

      return
      end


c ---------------------------------------------------------------------
      subroutine calc_coszen (jday)
c ---------------------------------------------------------------------
c
c This subroutine calculates coszen for every latitude and every sub
c -daily time step, storing the results in coszen_pre(lat, istep). It
c also calculates the daily mean of coszen for every latitude, storing
c the results in coszen_mean(lat). For times when the sun is below the
c horizon, coszen is set to 0.
c
c WJS 01.29.10: Some design notes:
c 
c The code here used to be inline in subroutine 'diurnal' in weather.f,
c through 01.29.10. It was moved here when adding the option to read
c radiation directly. In this case, we need to know the daily integral
c of coszen in order to convert from daily radiation to daily
c transmission fraction, so we need to precompute coszen for every time
c step in the day.
c
c For efficiency (in time and memory) we just do this calculation once
c per latitude band, rather than once per grid cell. This is fine as
c long as 'time' in the model refers to local time (e.g., time = 1 means
c 1 AM local time for all grid cells). However, this would need to be
c changed if 'time' referred to, say, Greenwich time, as might be the
c case for a run coupled to an atmoshperic model. In that case, coszen
c would depend on longitude as well as latitude, for a given time.
c However, for a coupled run like this, we should have direct access to
c the incoming radiation in this time step (rather than just daily
c radiation), which would remove the need for the integral of coszen. So
c we could move the code back into subroutine 'diurnal' in this case.
c (If we had incoming radiation in each time step, the code could be
c analogous to the single-site code, in diurnalmet.)
c
      use comatm, only : coszen_pre, coszen_mean
      use comgrid
      use compar
      use comwork
c
      implicit none
c
c Arguments
c
      integer jday    ! current day
c
c Local variables
c
      integer 
     >  lat,          ! looping variable over latitudes
     >  niter,        ! number of sub-daily time steps per day
     >  istep         ! looping variable over sub-daily time steps

      real 
     >  orbit,        ! earth's orbital angle (around the sun) in radians
     >  xdecl,        ! solar declination angle
     >  time,         ! time in seconds
     >  rtime,        ! time in hours
     >  angle,        ! solar hour angle in radians
     >  xlat          ! latitude in radians
c 
c Externals
c
      real calc_orbit, calc_xdecl, calc_angle
      real get_time

c coszen_mean will first store sums; initialize sums here
      do lat = 1, nlatsub
        coszen_mean(lat) = 0.0
      end do
      
c Calculate parameters that are constant in space and time for this day
      orbit = calc_orbit(jday)
      xdecl = calc_xdecl(orbit)

c Loop over time and latitude
      niter = int (86400.0 / dtime)
      do istep = 1, niter
        time = get_time(istep)
        rtime = time / 3600.0
        angle = calc_angle(rtime)

        do lat = 1, nlatsub
          xlat = latscale(lat) * pi / 180.0

c calculate the cosine of the solar zenith angle for the given latitude,
c day and time; store result in coszen_pre ('precomputed coszen')
          coszen_pre(lat, istep) = 
     >      max (0.0, (sin(xlat) * sin(xdecl) +
     >                 cos(xlat) * cos(xdecl) * cos(angle)))
        
          coszen_mean(lat) = coszen_mean(lat) + coszen_pre(lat, istep)  ! sums for now
        end do
      end do

c convert sums to means
      do lat = 1, nlatsub
        coszen_mean(lat) = coszen_mean(lat) / niter
      end do      

      return
      end



