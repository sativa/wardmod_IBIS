

c ---------------------------------------------------------------------
      module compar 
c     Last edited by Bill Sacks 03.24.10
c ---------------------------------------------------------------------
c
      use comgrid

      implicit none
      save

c ------------------------
c main parameters for ibis
c ------------------------
c WJS (02.11.09): grid-related parameters (nlon, nlat, npoi, xres, yres)
c moved to griddef.h
c        
      integer nband,   ! number of solar radiation wavebands
     >        nsoilay, ! number of soil layers
     >        nsnolay, ! number of snow layers
     >        npft,    ! number of total plant functional types, including crops
     >        scpft,   ! starting index for crop pfts - C. Kucharik 
     >        ecpft    ! ending index for crop pfts 
c
c WJS (03.25.10): Added c_to_k, but so far I'm only using it in new code, not in existing code
c 
      real pi,          ! you know, that constant thingy
     >     c_to_k,      ! conversion from degrees C to K
     >     un_init      ! value denoting uninitialized data
c
c -------------------------------
c state description configuration
c -------------------------------
c
      parameter (nband   = 2,
     >           nsoilay = 11,
     >           nsnolay = 3,
     >           scpft   = 13,
     >           ecpft   = 15,
     >           npft    = 15)
c
c ------------------------------------------------------------
c PFTs (so we can refer to them by name rather than by number)
c ------------------------------------------------------------
c WJS (03.24.10): I added these for later use, but so far they are unused in existing code

      integer, parameter :: 
     >  trop_bl_evergreen = 1,  ! tropical broadleaf evergreen trees
     >  trop_bl_ddecid = 2,     ! tropical broadleaf drought-deciduous trees
     >  wtemp_bl_evergreen = 3, ! warm-temperate broadleaf evergreen trees
     >  temp_con_evergreen = 4, ! temperate conifer evergreen trees
     >  temp_bl_cdecid = 5,     ! temperate broadleaf cold-deciduous trees
     >  bor_con_evergreen = 6,  ! boreal conifer evergreen trees
     >  bor_bl_cdecid = 7,      ! boreal broadleaf cold-deciduous trees
     >  bor_con_cdecid = 8,     ! boreal conifer cold-deciduous trees
     >  evergreen_shrub = 9,    ! evergreen shrubs
     >  cdecid_shrub = 10,      ! cold-deciduous shrubs
     >  c4_grass = 11,          ! warm (c4) grasses
     >  c3_grass = 12,          ! cool (c3) grasses
     >  soy = 13,               ! soybeans
     >  maize = 14,             ! maize
     >  wheat = 15              ! wheat

c --------------
c some constants
c --------------
c
      parameter (pi = 3.1415927)
      parameter (c_to_k = 273.16)
      parameter (un_init = huge(0.))
c
c ---------------------
c other misc. constants
c ---------------------
c
c tmax_wt and tmin_wt: For calculating tave from tmin and tmax, or tmin and tmax from tave and trange.
c Some relationships are (where trange = tmax - tmin):
c   tave = tmax_wt * tmax + tmin_wt * tmin
c   tmax = (tave - tmin_wt * tmin) / tmax_wt
c   tmax = tave + tmin_wt * trange
c   tmin = (tave - tmax_wt * tmax) / tmin_wt
c   tmin = tave - tmax_wt * trange
c Note that tmax_wt is equivalent to the daily mean value of gamma
      real
     >  tmax_wt,       ! weighting of tmax when calculating tave
     >  tmin_wt        ! weighting of tmin when calculating tave

      parameter (tmax_wt = 0.5,
     >           tmin_wt = 1 - tmax_wt)
c
c ------------------------------------------
c various constants, initialized in the code
c ------------------------------------------
c
      real 
     >  epsilon,       ! small quantity to avoid zero-divides and other
     >                 ! truncation or machine-limit troubles with small
     >                 ! values. should be slightly greater than o(1)
     >                 ! machine precision
     >  dtime,         ! model timestep (seconds)
     >  stef,          ! stefan-boltzmann constant (W m-2 K-4)
     >  vonk,          ! von karman constant (dimensionless)
     >  grav,          ! gravitational acceleration (m s-2)
     >  tmelt,         ! freezing point of water (K)
     >  hfus,          ! latent heat of fusion of water (J kg-1)
     >  hvap,          ! latent heat of vaporization of water (J kg-1)
     >  hsub,          ! latent heat of sublimation of ice (J kg-1)
     >  ch2o,          ! specific heat of liquid water (J deg-1 kg-1)
     >  cice,          ! specific heat of ice (J deg-1 kg-1)
     >  cair,          ! specific heat of dry air at constant pressure (J deg-1 kg-1)
     >  cvap,          ! specific heat of water vapor at constant pressure (J deg-1 kg-1)
     >  rair,          ! gas constant for dry air (J deg-1 kg-1)
     >  rvap,          ! gas constant for water vapor (J deg-1 kg-1)
     >  cappa,         ! rair/cair
     >  rhow           ! density of liquid water (all types) (kg m-3)
c
      real 
     >  garea(npoi),   ! area of each gridcell (m**2)
     >  vzero(npoi)    ! a real array of zeros, of length npoi
c
      integer 
     >  ndaypy,        ! number of days per year
     >  nlonsub,       ! number of longitude points for subsetting
     >  nlatsub        ! number of latitude points for subsetting
c
      integer
     >  ndaypm(12)     ! number of days per month
c
c
c ---------------------------------------------------------------
c Global variables for branch adamward-THMB (Y.Li 2014-06-23)
c ---------------------------------------------------------------
c
      integer
     > isoil_layer,    ! 0: no; 1: yes    Read-in soils with layered data rom a NetCDF file
     > iyr_pft,        ! 0: no; 1: yes    Read-in land use on a ear-by-year basis from NetCDF files
     > iyr_fert,       ! 0: no; 1: yes    Read-in fertilizer application n a year-by-year basis from NetCDF files
     > ihourlyout,     ! 0: no hourly output, 1: hourly output
     > iagro_thmb      ! 0: no coupled agro-ibis thmb simulation(default); 1: yes

      character*300::   agro_thmb_path  !path for communication files for coupled agro-ibis thmb

      integer:: it_hour   !current hour of a day
      integer:: commu_MaxWait, commu_sleep
      integer:: first_day  !1:yes; 0:no  first day of simulation?
      integer:: icday      !cummulative days of a year

      real,parameter:: stamp_eps=0.001
      real
     > thmb_lakem(npoi)  !input from thmb, assgin to variable <zwpud> at <soil.f> at first hour of a day

      !--file IO
      integer,parameter:: thmb_FnameID=1001,agro_FnameID=1002
      character*40,parameter:: thmb_Fname ='thmb.communicate',
     >                         agro_Fname = 'agro.communicate'

      !--error handler
      logical:: checked_grid !check grid consistency for THMB and Agro-IBIS; check only once
      logical:: error_showed_sumday !show only once for error   in 'ERROR in sumday'
      logical:: error_showed_plen   !show only once for warning in 'WARNING: plen changed'
      real,parameter:: grid_torlerance=0.01  !grid tolerance for THMB and Agro-IBIS

      !--variables for user-defined grid mode
      integer:: flg_udg
      integer:: udg_is,udg_ie,udg_js,udg_je
      real(kind=4):: lon_llcorner,lat_llcorner

      !--other
      integer:: flg_debug,debug_ncounts
      integer:: flg_detail
      integer:: flg_daily4thmb
      integer:: flg_sandlayers
      integer:: flg_dailyin_freq
      integer:: flg_dailyout_freq
      integer:: flg_cmask
      integer:: flg_coverCrops
      integer:: flg_mngt
      integer:: dyn_idx  !counter for dailyout_freq

      integer:: irestart_glo   !same as irestart, used just to avoid pass the argument
      integer:: irstyear_glo   !same as iyear0,   used just to avoid pass the argument
      integer:: iyear0_glo     !same as iyear0,   used just to avoid pass the argument
      integer:: nrun_glo       !same as iyear0,   used just to avoid pass the argument
      integer:: imonthout_glo  !same as imonthout
      integer:: iyearout_glo   !same as iyearout
      integer:: iyrlast_glo    !same as iyrlast
      integer:: idailyout_glo  !same as idailyout_glo


      character*200 fileName_cmask,varName_cmask
      character*200 myprocDir  !base directory for proc<n>
      integer:: flg_mpi        !run in mpi? 1=yes, 0=no (default)
      integer:: istart_all,iend_all,jstart_all,jend_all
      integer:: flg_ghost_io   !ghost IO? 1=yes, 0=no (default)
      ! S. Slavin added flg_rdopt. Set in main.f and used by ies-io.f
      integer :: flg_rdopt = 0  !flag for setting read optimization on/off

      !--nc file IO
      integer:: nc_file_bufrsize_KB        !buffer size (KB) for nc file IO, default = 8 
      integer,parameter:: size_1KB = 1024  !1 KB

      end module compar
