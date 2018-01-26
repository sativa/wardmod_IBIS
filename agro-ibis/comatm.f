c
c ---------------------------------------------------------------------
      module comatm
c ---------------------------------------------------------------------
c
c In addition to data, this module includes a subroutine to allocate
c space for all allocatable variables in the module.
c
      use comgrid
      use compar, only : nband

      implicit none
      save
c      
      real 
     >  coszen(npoi),      ! cosine of solar zenith angle
     >  fira(npoi)         ! incoming ir flux (W m-2)
c
      real 
     >  solad(npoi,nband), ! direct downward solar flux (W m-2)
     >  solai(npoi,nband), ! diffuse downward solar flux (W m-2)
     >  asurd(npoi,nband), ! direct albedo of surface system
     >  asuri(npoi,nband)  ! diffuse albedo of surface system 
c
      real 
     >  ua(npoi),          ! wind speed (m s-1)
     >  ta(npoi),          ! air temperature (K)
     >  qa(npoi),          ! specific humidity (kg_h2o/kg_air)
     >  raina(npoi),       ! rainfall rate (mm/s or kg m-2 s-1)
     >  rh(npoi),          ! relative humidity(%)
     >  snowa(npoi)        ! snowfall rate (mm/s or kg m-2 s-1 of water)
c
      real 
     >  psurf(npoi),       ! surface pressure (Pa)
     >  cloud(npoi),       ! cloud fraction
     >  rads(npoi),        ! solar radiation at the surface (W m-2)
     >  trans(npoi),       ! solar transmission through the atmosphere (fraction)
     >  td(npoi),          ! daily average temperature (K)
     >  tmax(npoi),        ! maximum daily temperature (K)
     >  tmin(npoi),        ! minimum daily temperature (K)
     >  qd(npoi),          ! daily average specific humidity (kg_h2o/kg_air)
     >  ud(npoi),          ! daily average wind speed (m/sec)
     >  precip(npoi),      ! daily precitation (mm/day)
     >  precipday(npoi,31),       
     >  precipdaysum(npoi)       
c
      real 
     >  xstore(npoi,3)     ! weather generator 'memory' matrix
c
      integer 
     >  iwet(npoi),        ! wet day / dry day flag
     >  iwetday(npoi,31), 
     >  iwetdaysum(npoi) 
c
      real 
     >  co2conc,           ! co2 concentration (mol/mol)
     >  o2conc             ! o2 concentration (mol/mol)
c
c--------------------------------------------------------------------------------
c added for Green-Ampt
      real
     >    t_sec,     ! time since midnight in second
     >    t_startp   ! time since midnight in second when precipitation starts 
c
c -------------------------------------------------------------------------------

c coszen_pre(lat, istep): Pre-computed values of coszen for every
c latitude and every sub-daily time step for the given day.
      real, dimension(:,:), allocatable :: coszen_pre ! size: nlatsub x niter

c coszen_mean(lat): Daily mean of coszen for every latitude for the given day.
      real, dimension(:), allocatable :: coszen_mean   ! size: nlatsub



      contains

c ---------------------------------------------------------------------
      subroutine comatm_init(niter)
c ---------------------------------------------------------------------
c
c Allocate dynamically-allocated arrays in comatm.
c
c This must be called after nlatsub is set.
c
      use comwork
c
c Arguments
c 
      integer niter   ! number of time iterations per day

      allocate(coszen_pre(nlatsub, niter))
      allocate(coszen_mean(nlatsub))      
      
      return
      end subroutine comatm_init
      

      end module comatm
