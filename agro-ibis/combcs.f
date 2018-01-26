c
c ---------------------------------------------------------------------
      module combcs
c ---------------------------------------------------------------------
c
      use comgrid

      implicit none
      save

      real
     >  xintopo(npoi),    ! topography (m)
     >  xinveg(npoi),     ! fixed vegetation map
     >  deltat(npoi)      ! absolute minimum temperature -
     >                    ! temp on average of coldest month (C)
c
      integer
     >  lmask(nlon,nlat)  ! landmask 0=water, 1=land

      real
     >  cmask(npoi)         !computational mask; 1 = cell to compute; 0 = no computation

c
      real
     >  xint(npoi,12),    ! climatological temp + anomaly (C)
     >                    !  Note: xint is either read in or set based on xintmin and xintmax.
     >  xintmin(npoi,12), ! climatoligical daily min temp. + anomaly (C)
     >                    !  Note: xintmin is either read in or set based on xint and xintrng.
     >  xintmax(npoi,12), ! climatoligical daily max temp. + anomaly (C)
     >                    !  Note: xintmax is either read in or set based on xint and xintrng.
     >  xinq(npoi,12),    ! climatological relative humidity + anomaly (%)
     >  xinprec(npoi,12), ! climatological precipition + anomaly (mm/day)
     >  xinwind(npoi,12), ! climatological wind speed + anomaly (m s-1)
     >  xincld(npoi,12),  ! climatological cloudiness + anomaly(%)
     >                    !  Note: xincld is only available if read_radiation is FALSE;
     >                    !  if read_radiation is TRUE, all values are set to un_init.
     >  xinrads(npoi,12), ! climatological solar radiation + anomaly (W m-2)
     >                    !  Note: xinrads is only available if read_radiation is TRUE;
     >                    !  if read_radiation is FALSE, all values are set to un_init.
     >  xinwet(npoi,12),  ! climatological wet days + anomaly (days/month)
     >                    !  Note: xinwet can be missing, in which case all values are set to un_init.
     >                    !  In this case we can't use the weather generator.
     >  xintrng(npoi,12)  ! climatological temp range + anomaly(C)
                          !  Note: xintrng is either read in or set based on xintmin and xintmax.
c
      real
     >  clmt(npoi,12),    ! climatological temperature (C)
     >                    !  Note: clmt is either read in or set based on clmtmin and clmtmax.
     >  clmtmin(npoi,12), ! climatological daily min temp. (C)
     >                    !  Note: clmtmin is either read in or set based on clmt and clmtrng.
     >  clmtmax(npoi,12), ! climatological daily min temp. (C)
     >                    !  Note: clmtmax is either read in or set based on clmt and clmtrng.
     >  clmq(npoi,12),    ! climatological relative humidity (%)
     >  clmprec(npoi,12), ! climatological precipitation (mm/day)
     >  clmw(npoi,12),    ! climatological wind speed (m s-1)
     >  clmwet(npoi,12),  ! climatological wet days (days/month)
     >                    !  Note: clmwet can be missing, in which case all values are set to un_init.
     >                    !  (see note on xinwet, above)
     >  clmcld(npoi,12),  ! climatological cloudiness (%)
     >                    !  Note: clmcld is only available if read_radiation is FALSE;
     >                    !  if read_radiation is TRUE, all values are set to un_init.
     >  clmrads(npoi,12), ! climatological solar radiation (W m-2)
     >                    !  Note: clmrads is only available if read_radiation is TRUE;
     >                    !  if read_radiation is FALSE, all values are set to un_init.
     >  clmtrng(npoi,12)  ! climatological temp range (C)
                          !  Note: clmtrng is either read in or set based on clmtmin and clmtmax.
c
c
c WJS (01.19.10): For the climate variables given below (not xinavepdate
c or xinavehybrid), these variables represent the input values from the
c daily climate files. This can either mean daily anomalies (if we are
c using a climate data set for which we apply daily anomalies to monthly
c values), or it can mean the daily values themselves (if we are using a
c climate data set that specifies the daily values directly).
      real
     >  xintd(npoi),      ! daily temperature anomaly (C, additive anomaly)
     >                    !  OR daily temperature (C)
     >                    !  Note: xintd is only available if read_tmin_tmax is FALSE;
     >                    !  if read_tmin_tmax is TRUE, all values are set to un_init
     >                    !  (in contrast to xint, xintd is NOT set based on xintmind and xintmaxd).
     >  xinqd(npoi),      ! daily relative humidity anomaly (fraction, multiplicative anomaly)
     >                    !  OR daily relative humidity (%)
     >  xinprecd(npoi),   ! daily precipitation anomaly (fraction, multiplicative anomaly)
     >                    !  OR daily precipitation (mm/day)
     >  xinwindd(npoi),   ! daily wind speed anomaly (fraction, multiplicative anomaly)
     >                    !  OR daily wind speed (m/s)
     >  xincldd(npoi),    ! daily cloud fraction anomaly (%, additive anomaly)
     >                    !  OR daily cloud fraction (%)
     >                    !  Note: xincldd is only available if read_radiation is FALSE;
     >                    !  if read_radiation is TRUE, all values are set to un_init.
     >  xinradsd(npoi),   ! daily solar radiation anomaly (units to be determined)
     >                    !  OR daily solar radiation (W m-2)
     >                    !  Note: xinradsd is only available if read_radiation is TRUE;
     >                    !  if read_radiation is FALSE, all values are set to un_init.
     >  xintrngd(npoi),   ! daily temp range anomaly (C, multiplicative anomaly)
     >                    !  OR daily temp range (C)
     >                    !  Note: xintrngd is only available if read_tmin_tmax is FALSE;
     >                    !  if read_tmin_tmax is TRUE, all values are set to un_init
     >                    !  (in contrast to xintrng, xintrngd is NOT set based on xintmind and xintmaxd).
     >  xintmaxd(npoi),   ! daily temp maximum anomaly (C, additive anomaly)  
     >                    !  OR daily temp maximum (C)
     >                    !  Note: xintmaxd is only available if read_tmin_tmax is TRUE;
     >                    !  if read_tmin_tmax is FALSE, all values are set to un_init
     >                    !  (in contrast to xintmax, xintmaxd is NOT set based on xintd and xintrngd).
     >  xintmind(npoi),   ! daily 2m temp minimum anomaly (C, additive anomaly)  
     >                    !  OR daily temp minimum (C)
     >                    !  Note: xintind is only available if read_tmin_tmax is TRUE;
     >                    !  if read_tmin_tmax is FALSE, all values are set to un_init
     >                    !  (in contrast to xintmin, xintmind is NOT set based on xintd and xintrngd).
     >  xinavepdate(npoi),! average crop planting date input
     >  xinavehybrid(npoi)! average corn crop gdd hybrid for period
c
      end module combcs
