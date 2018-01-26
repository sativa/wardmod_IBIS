!
!***********************************************************************
!
! - Dynamic IO
!
!                       Yuwei Li 2015-12-21
!  
!***************************** Note ************************************
!  
module daily_variables

  implicit NONE

  save

  !parameters
  integer,parameter:: max_days_per_month = 31, max_days_per_year = 366

  !variables
  !dimension(npoi,max_days_per_*)
  real(kind=4),dimension(:,:),allocatable::   adrain_dyn,adtrunoff_dyn,adsrunoff_dyn,addrainage_dyn,adwsoi_dyn,adwisoi_dyn,adsnod_dyn, &
                                              adsnof_dyn,adco2ratio_dyn,adco2mic_dyn,adglail_dyn,adblail_dyn,dnileach_dyn,aplantn_dyn,ftot_dyn, &
                                              drntot_dyn,adtrans_dyn,adevap_dyn
  !dimension(npoi,nsoilay,max_days_per_*)
  real(kind=4),dimension(:,:,:),allocatable:: csoln_dyn

  !dimension(npoi,npft,max_days_per_*)
  real(kind=4),dimension(:,:,:),allocatable:: plai_dyn,stressn_dyn,totnuptake_dyn,fnleaf_dyn,fnstem_dyn,fnroot_dyn,fngrain_dyn,fnplant_dyn,tnplant_dyn

  !dimension(npoi,2,max_days_per_*)
  real(kind=4),dimension(:,:,:),allocatable:: lai_dyn

  !**climate input variables
  !dimension(npoi,max_days_per_*)
  real(kind=4),dimension(:,:),allocatable:: xinprecd_dyn,xintmind_dyn,xintmaxd_dyn,xintd_dyn,xintrngd_dyn,xinradsd_dyn, &
                                            xincldd_dyn ,xinwindd_dyn,xinqd_dyn

  !**daily4thmb variables
  !dimension(npoi,nsoilay,max_days_per_*) -- daily-scale
  real(kind=4),dimension(:,:,:),allocatable::  addrainage_layer_dyn,adtotnleach_layer_dyn

  !dimension(npoi,24,max_days_per_*) -- hourly-scale
  real(kind=4),dimension(:,:,:),allocatable::   adsrunoff_hourly_dyn,adrain_hourly_dyn,adevap_hourly_dyn
  !dimension(npoi,nsoilay,24,max_days_per_*)
  real(kind=4),dimension(:,:,:,:),allocatable:: addrainage_layer_hourly_dyn,adtotnleach_layer_hourly_dyn



end module daily_variables
!
!***********************************************************************
!
! - Management modules/functionalities
!
!                       Yuwei Li 2016-04-11
!  
!***************************** Note ************************************
! 
module management_para

  implicit NONE

  save

  !<management.nml>
  !namelist-climate
  character*10 climate_mngt_mode
  character*300 climate_mngt_file

  !namelist-fertilizer
  character*10 fertilizer_mngt_mode
  integer:: fertilizer_mngt_year_start,fertilizer_mngt_year_end
  character*300 fertilizer_mngt_file

  integer,parameter:: MAX_YEARS = 200

  !climate variables
  integer:: climateMngt_yearsList(MAX_YEARS)

  !fertilizer variables
  integer:: fertMngt_ys_idx,fertMngt_ye_idx,fertMngt_nyears

  !dimension(npoi) (fertilizer_mngt_mode='dat') or (1) otherwise
  real(kind=4),dimension(:),allocatable::  mngt_fertnitro0_13,mngt_fertnitro0_14,mngt_fertnitro0_15  !total fertilizer in a year;used to distribute to daily per management

  !dimension(366+1) (fertilizer_mngt_mode='dat') or (1) otherwise; index 1 is for the date at planting;
  real(kind=4),dimension(:),allocatable:: fertMngt_pct  !range = 0~1, daily portion of fertilizer to apply out of the one-year fertilizer

end module management_para
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : alloc_dailyout_variables
! purpose : allocate and initialize variables for daily output
!
! note: called only when flg_dailyout_freq != 0
!
!                           Y.Li 2015-12-21
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine alloc_dailyout_variables()
  use compar, only: flg_dailyout_freq,npft,nsoilay
  use comgrid,only: npoi
  use daily_variables

  implicit NONE
  !--local
  integer ierr,dyn_dim

  write(6,'(A)',advance='no') 'allocating dailyout variables...'
 
  if(flg_dailyout_freq == 1) then      !one month
    dyn_dim = max_days_per_month
  else if(flg_dailyout_freq == 2) then !one year
    dyn_dim = max_days_per_year  
  end if

  !dimension(npoi,max_days_per_*)
  allocate(adrain_dyn    (npoi,dyn_dim),STAT=ierr)
  allocate(adtrunoff_dyn (npoi,dyn_dim),STAT=ierr)
  allocate(adsrunoff_dyn (npoi,dyn_dim),STAT=ierr)
  allocate(addrainage_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(adwsoi_dyn    (npoi,dyn_dim),STAT=ierr)
  allocate(adwisoi_dyn   (npoi,dyn_dim),STAT=ierr)
  allocate(adsnod_dyn    (npoi,dyn_dim),STAT=ierr)

  allocate(adsnof_dyn    (npoi,dyn_dim),STAT=ierr)
  allocate(adco2ratio_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(adco2mic_dyn  (npoi,dyn_dim),STAT=ierr)
  allocate(adglail_dyn   (npoi,dyn_dim),STAT=ierr)  !green lai of lower canopy
  allocate(adblail_dyn   (npoi,dyn_dim),STAT=ierr)  !brown lai of lower canopy
  allocate(dnileach_dyn  (npoi,dyn_dim),STAT=ierr)
  allocate(aplantn_dyn   (npoi,dyn_dim),STAT=ierr)
  allocate(ftot_dyn      (npoi,dyn_dim),STAT=ierr)
  allocate(drntot_dyn    (npoi,dyn_dim),STAT=ierr)

  allocate(adtrans_dyn   (npoi,dyn_dim),STAT=ierr)
  allocate(adevap_dyn    (npoi,dyn_dim),STAT=ierr)

  !dimension(npoi,nsoilay,max_days_per_*)
  allocate(csoln_dyn            (npoi,nsoilay,dyn_dim),STAT=ierr)

  !dimension(npoi,npft,max_days_per_*)
  allocate(plai_dyn      (npoi,npft,dyn_dim),STAT=ierr)

  allocate(stressn_dyn   (npoi,npft,dyn_dim),STAT=ierr)
  allocate(totnuptake_dyn(npoi,npft,dyn_dim),STAT=ierr)
  allocate(fnleaf_dyn    (npoi,npft,dyn_dim),STAT=ierr)
  allocate(fnstem_dyn    (npoi,npft,dyn_dim),STAT=ierr)
  allocate(fnroot_dyn    (npoi,npft,dyn_dim),STAT=ierr)
  allocate(fngrain_dyn   (npoi,npft,dyn_dim),STAT=ierr)
  allocate(fnplant_dyn   (npoi,npft,dyn_dim),STAT=ierr)
  allocate(tnplant_dyn   (npoi,npft,dyn_dim),STAT=ierr)

  !dimension(npoi,2,max_days_per_*)
  allocate(lai_dyn (npoi,2,dyn_dim),STAT=ierr)

  write(6,'(A)') 'done'

  return
end subroutine alloc_dailyout_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : alloc_dailyin_variables
! purpose : allocate and initialize variables for daily input
!
! note: called only when flg_dailyin_freq != 0
!
!                           Y.Li 2016-01-19
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine alloc_dailyin_variables()
  use compar, only: flg_dailyin_freq
  use comgrid,only: npoi
  use daily_variables

  implicit NONE
  !--local
  integer ierr,dyn_dim

  write(6,'(A)',advance='no') 'allocating dailyin variables...'

  if(flg_dailyin_freq == 1) then      !one month
    dyn_dim = max_days_per_month
  else if(flg_dailyin_freq == 2) then !one year
    dyn_dim = max_days_per_year  
  end if

  !**climate input variables
  !dimension(npoi,max_days_per_*)
  allocate(xinprecd_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(xintmind_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(xintmaxd_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(xintd_dyn   (npoi,dyn_dim),STAT=ierr)
  allocate(xintrngd_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(xinradsd_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(xincldd_dyn (npoi,dyn_dim),STAT=ierr)
  allocate(xinwindd_dyn(npoi,dyn_dim),STAT=ierr)
  allocate(xinqd_dyn   (npoi,dyn_dim),STAT=ierr)

  write(6,'(A)') 'done'

  return
end subroutine alloc_dailyin_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : alloc_daily4thmb_variables
! purpose : allocate and initialize variables specifically for daily4thmb output
!
! note: called only when flg_dailyout_freq != 0 and flg_daily4thmb = 1 or 2
!
!                           Y.Li 2016-01-22
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine alloc_daily4thmb_variables()
  use compar, only: flg_dailyout_freq,flg_daily4thmb,nsoilay,dtime
  use comgrid,only: npoi
  use daily_variables

  implicit NONE
  !--local
  integer ierr,dyn_dim

  write(6,'(A)',advance='no') 'allocating daily4thmb variables...'
 
  if(flg_dailyout_freq == 1) then      !one month
    dyn_dim = max_days_per_month
  else if(flg_dailyout_freq == 2) then !one year
    dyn_dim = max_days_per_year  
  end if

  if(flg_daily4thmb == 1) then !daily-scale

    !dimension(npoi,nsoilay,max_days_per_*)
    allocate(addrainage_layer_dyn  (npoi,nsoilay,dyn_dim),STAT=ierr)
    allocate(adtotnleach_layer_dyn (npoi,nsoilay,dyn_dim),STAT=ierr)

  elseif(flg_daily4thmb == 2) then  !hourly-scale

    !simple check
    if( int(86400.0/dtime) /= 24 ) then !not 24 hours (steps) per day
      write(6,*) '(alloc_daily4thmb_variables)Error! expecting 24 steps per day,i.e. dtime = 3600'
      write(6,*) 'Stopping the program...'
      stop 11
    end if

    !dimension(npoi,24,max_days_per_*) 
    allocate(adsrunoff_hourly_dyn (npoi,24,dyn_dim),STAT=ierr)
    allocate(adrain_hourly_dyn    (npoi,24,dyn_dim),STAT=ierr)
    allocate(adevap_hourly_dyn    (npoi,24,dyn_dim),STAT=ierr)

    !dimension(npoi,nsoilay,24,max_days_per_*) 
    allocate(addrainage_layer_hourly_dyn  (npoi,nsoilay,24,dyn_dim),STAT=ierr)
    allocate(adtotnleach_layer_hourly_dyn (npoi,nsoilay,24,dyn_dim),STAT=ierr)

  end if  !on if flg_daily4thmb
  
  write(6,'(A)') 'done'

  return
end subroutine alloc_daily4thmb_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : dealloc_dailyout_variables
! purpose : deallocate variables for daily output
!
! note: called only when flg_dailyout_freq != 0
!
!                           Y.Li 2015-12-21
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine dealloc_dailyout_variables()
  use daily_variables
  implicit NONE
  !--local
  integer ierr

  !dimension(npoi,max_days_per_*)
  deallocate(adrain_dyn    ,STAT=ierr)
  deallocate(adtrunoff_dyn ,STAT=ierr)
  deallocate(adsrunoff_dyn ,STAT=ierr)
  deallocate(addrainage_dyn,STAT=ierr)
  deallocate(adwsoi_dyn    ,STAT=ierr)
  deallocate(adwisoi_dyn   ,STAT=ierr)
  deallocate(adsnod_dyn    ,STAT=ierr)

  deallocate(adsnof_dyn    ,STAT=ierr)
  deallocate(adco2ratio_dyn,STAT=ierr)
  deallocate(adco2mic_dyn  ,STAT=ierr)
  deallocate(adglail_dyn   ,STAT=ierr)
  deallocate(adblail_dyn   ,STAT=ierr)
  deallocate(dnileach_dyn  ,STAT=ierr)
  deallocate(aplantn_dyn   ,STAT=ierr)
  deallocate(ftot_dyn      ,STAT=ierr)
  deallocate(drntot_dyn    ,STAT=ierr)

  deallocate(adtrans_dyn   ,STAT=ierr)
  deallocate(adevap_dyn    ,STAT=ierr)

  !dimension(npoi,nsoilay,max_days_per_*)
  deallocate(csoln_dyn            ,STAT=ierr)

  !dimension(npoi,npft,max_days_per_*)
  deallocate(plai_dyn      ,STAT=ierr)
  deallocate(stressn_dyn   ,STAT=ierr)
  deallocate(totnuptake_dyn,STAT=ierr)
  deallocate(fnleaf_dyn    ,STAT=ierr)
  deallocate(fnstem_dyn    ,STAT=ierr)
  deallocate(fnroot_dyn    ,STAT=ierr)
  deallocate(fngrain_dyn   ,STAT=ierr)
  deallocate(fnplant_dyn   ,STAT=ierr)
  deallocate(tnplant_dyn   ,STAT=ierr)

  !dimension(npoi,2,max_days_per_*)
  deallocate(lai_dyn ,STAT=ierr)

  return
end subroutine dealloc_dailyout_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : dealloc_dailyin_variables
! purpose : deallocate variables for daily input
!
! note: called only when flg_dailyin_freq != 0
!
!                           Y.Li 2016-01-19
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine dealloc_dailyin_variables()
  use daily_variables
  implicit NONE
  !--local
  integer ierr

  !**climate input variables
  !dimension(npoi,max_days_per_*)
  deallocate(xinprecd_dyn,STAT=ierr)
  deallocate(xintmind_dyn,STAT=ierr)
  deallocate(xintmaxd_dyn,STAT=ierr)
  deallocate(xintd_dyn   ,STAT=ierr)
  deallocate(xintrngd_dyn,STAT=ierr)
  deallocate(xinradsd_dyn,STAT=ierr)
  deallocate(xincldd_dyn ,STAT=ierr)
  deallocate(xinwindd_dyn,STAT=ierr)
  deallocate(xinqd_dyn   ,STAT=ierr)

  return
end subroutine dealloc_dailyin_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : dealloc_daily4thmb_variables
! purpose : deallocate variables for daily4thmb output
!
! note: called only when flg_dailyout_freq != 0 and flg_daily4thmb = 1 or 2
!
!                           Y.Li 2016-01-22
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine dealloc_daily4thmb_variables()
  use compar, only: flg_daily4thmb
  use daily_variables
  implicit NONE
  !--local
  integer ierr

  if(flg_daily4thmb == 1) then
    !dimension(npoi,nsoilay,max_days_per_*)
    deallocate(addrainage_layer_dyn  ,STAT=ierr)
    deallocate(adtotnleach_layer_dyn ,STAT=ierr)

  elseif(flg_daily4thmb == 2) then

    !dimension(npoi,24,max_days_per_*) 
    deallocate(adsrunoff_hourly_dyn ,STAT=ierr)
    deallocate(adrain_hourly_dyn    ,STAT=ierr)
    deallocate(adevap_hourly_dyn    ,STAT=ierr)

    !dimension(npoi,nsoilay,24,max_days_per_*) 
    deallocate(addrainage_layer_hourly_dyn  ,STAT=ierr)
    deallocate(adtotnleach_layer_hourly_dyn ,STAT=ierr)

  end if  !on if flg_daily4thmb

  return
end subroutine dealloc_daily4thmb_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : write_dailyout_dyn
! purpose : write dailyout at the end of each month or year
!
! note:
!    1) called only when flg_dailyout_freq != 0 and at the end of the frquency period
!    2) output daily or hourly for thmb when flg_daily4thmb = 1 or 2
!        
!
!                           Y.Li 2015-12-21
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine write_dailyout_dyn(imonth,iyear)

  use comwork, only: lonscale,latscale,OCEAN
  use comgrid, only: npoi
  use compar
  use comcrop,only: isoilay
  use daily_variables
  use netcdf_wrappers

  implicit NONE

  include 'netcdf.inc'

  integer:: iyear,imonth
  integer:: vardim(3),idx,istart_time,ndays_pastMonths
  integer:: ndims,n3rd,istart(4), icount(4) ! for writing vars
  integer:: idies,istat
  integer:: i,j,ilayer
  integer:: layer_s,layer_e

  character*10  cdate        ! date to use in history attribute in files
  character*10  tdate        ! character date for time step
  character*21  tunits       ! time units
  character*80  dimnames(4)  ! names of dimensions for vars
  character*200 filen        ! file name
  character*200 subfile      ! sub-file name, for some split file due to file size
  character*200 name3rd,long3rd,units3rd,pos3rd
  character*200 varName,longname
  character*80  units

  character*80  pftdef(npft)  ! plant functional type defs (not currently used)
  character*13  canopies(2)   ! canopy definitions

  character*2   num_lay       !char for soil layer

  real:: ftime(max_days_per_year),&              ! real form of nday 
         tweight(max_days_per_year)              ! number of days in daily average = 1. dummy (1)

  real:: vals3rd(100)         !large enough to serve possible nlevels (npfts,nsoilay, etc.)
  logical:: new_file          !Each file for each year
  logical:: is_hourly         !is daily4thmb hourly?
  character*4 num,num2

  !--initialization
  new_file = .true.
  if(flg_dailyout_freq == 1 .and. imonth /= 1) then
    new_file = .false.
  end if

  do i =1, 100
    vals3rd(i) = real(i)
  end do 

  !get cummulative days before current month
  ndays_pastMonths = 0
  if(imonth > 1) then
    do i = 1, imonth-1
      ndays_pastMonths = ndays_pastMonths + ndaypm(i)
    end do
  end if

  if(flg_dailyout_freq == 1) then      !one month
    istart_time = ndays_pastMonths + 1
  else if(flg_dailyout_freq == 2) then !one year
    istart_time = 1
  end if

  write(num,'(i4.4)')  iyear
  write(num2,'(i4.4)') iyear-1
  tunits    = 'days since '//num2//'-12-31 00:00:00'

  istart  = (/1,       1,       1, istart_time /)
  icount  = (/nlonsub, nlatsub, 1, dyn_idx     /)  !3rd icount may need to change depending on variable

  idx = 0
  do i = istart_time,istart_time+dyn_idx-1
    idx = idx + 1
    ftime(idx)   = real(i)
    tweight(idx) = 365.0
  end do

  ndims       =  4
  dimnames(1) = 'longitude'
  dimnames(2) = 'latitude'
  dimnames(3) = 'level'
  dimnames(4) = 'time'

  name3rd  = trim(dimnames(3))
  long3rd  = trim(name3rd)
  units3rd = ''
  pos3rd   = ''
  n3rd     = 1            !size of 3rd dimension (nlevel)

  vardim(1)   = npoi      !vardim could be 2d or 3d for the variable, here dimension of 3 is used for generality
  vardim(2)   = 1
  vardim(3)   = dyn_idx   !time dimension

  !<define char constants>
  pftdef(1)   = 'trbrevtr - tropical broadleaf evergreen trees'//char(0)
  pftdef(2)   = 'trbrdetr - tropical broadleaf drought-deciduous trees'//char(0)
  pftdef(3)   = 'wtbrevtr - warm-temperate broadleaf evergreen trees'//char(0)
  pftdef(4)   = 'tecoevtr - temperate conifer evergreen trees'//char(0)
  pftdef(5)   = 'tebrdetr - temperate broadleaf cold-deciduous trees'//char(0)
  pftdef(6)   = 'bocoevtr - boreal conifer evergreen trees'//char(0)
  pftdef(7)   = 'bocodetr - boreal conifer cold-deciduous trees'//char(0)
  pftdef(8)   = 'bobrdetr - boreal broadleaf cold-deciduous trees'//char(0)
  pftdef(9)   = 'evsh - evergreen shrubs'//char(0)
  pftdef(10)  = 'desh - deciduous shrubs'//char(0)
  pftdef(11)  = 'c4gr - warm (c4) grasses'//char(0)
  pftdef(12)  = 'c3gr - cool (c3) grasses'//char(0)
  pftdef(13)  = 'c3 crop - soybean'//char(0)
  pftdef(14)  = 'c4 crop - corn'//char(0)
  pftdef(15)  = 'c3 crop - wheat'//char(0)

  canopies(1) = 'lower canopy'//char(0)
  canopies(2) = 'upper canopy'//char(0)
  !</define char constants>

  !--begin the work

  !NOTE (Li): To output daily4thmb data only (but not other daily files), set idailyout (in ibis.infile) to 2
  !           A better design/structure is recommended if this feature is a common need 
  !           To erase this feature, get rid of codes related to idailyout_glo == 2 (idailyout_glo is global var for idailyout)

  if(idailyout_glo == 2) then  !skip all other daily output
    goto  1131
  end if
  
  !**rainfall (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'rain.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average rainfall','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'rain','average rainfall','mm/day',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'rain',vardim,istart,icount,adrain_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, rain'
     stop 1
  end if

  !**trunoff (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'trunoff.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average total runoff','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'trunoff','average total runoff','mm/day',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'trunoff',vardim,istart,icount,adtrunoff_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, trunoff'
     stop 1
  end if

  !**srunoff (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'srunoff.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average surface runoff','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'srunoff','average surface runoff','mm/day',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'srunoff',vardim,istart,icount,adsrunoff_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, srunoff'
     stop 1
  end if

  !**drainage (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'drainage.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average drainage','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'drainage','average drainage','mm/day',ndims,dimnames,OCEAN,istat)

    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'drainage',vardim,istart,icount,addrainage_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, drainage'
     stop 1
  end if

  !**wsoi (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'wsoi.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average soil moisture','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'wsoi','average soil moisture','fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'wsoi',vardim,istart,icount,adwsoi_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, wsoi'
     stop 1
  end if

  !**wisoi (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'wisoi.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average soil ice','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'wisoi','average soil ice','fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'wisoi',vardim,istart,icount,adwisoi_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, wisoi'
     stop 1
  end if

  !**snod (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'snod.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average snow depth','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'snod','average snow depth','meters',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'snod',vardim,istart,icount,adsnod_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, snod'
     stop 1
  end if

  !**snof (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'snof.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average snow fraction','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'snof','average snow fraction','fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'snof',vardim,istart,icount,adsnof_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, snof'
     stop 1
  end if

  !**co2ratio (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'co2ratio.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily average ratio of root to total soil co2 flux','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'co2ratio','average co2 ratio','fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'co2ratio',vardim,istart,icount,adco2ratio_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, co2ratio'
     stop 1
  end if

  !**co2mic (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'co2mic.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily flux of carbon due to soil microbe co2 flux','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'co2mic','soil microbe carbon flux','kg/m^2',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'co2mic',vardim,istart,icount,adco2mic_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, co2mic'
     stop 1
  end if

  !**lower canopy daily green & brown lai (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'lailower.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily green & brown lai-leaf area index of lower vegetation canopy','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'glail','daily green lai of lower canopy','m2/m2',ndims,dimnames,OCEAN,istat)
    call inivar(idies,'blail','daily brown lai of lower canopy','m2/m2',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'glail',vardim,istart,icount,adglail_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, glail'
     stop 1
  end if

  call write4dvar_ver2(filen,'blail',vardim,istart,icount,adblail_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, blail'
     stop 1
  end if

  !**daily rate of nitrate leaching (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'leachr.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'rate of nitrogen leached from profile (kg N m-2 y-1)','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'dnileach','rate of nitrogen leaching','kg n ha-1 y-1',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'dnileach',vardim,istart,icount,dnileach_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, dnileach'
     stop 1
  end if

  !**instantaneous plant available inorganic nitrogen  (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'aplantn.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'inorganic plant available inorganic nitrogen (kg N m-2 )','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'aplantn','inorganic plant available nitrogen','kg n m-2',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'aplantn',vardim,istart,icount,aplantn_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, aplantn'
     stop 1
  end if

  !**cumulative nitrogen leaching (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'cumnleach.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'accumulated nitrogen leached from profile at specific depth (kg N m-2)','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'ftot','cumulative nitrogen leaching','kg n ha-1',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'ftot',vardim,istart,icount,ftot_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, ftot'
     stop 1
  end if

  !**cumulative drainage (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'cumdrainage.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'accumulated drainage from profile at specific depth (mm h20)','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'drntot','cumulative drainage','mm h20',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'drntot',vardim,istart,icount,drntot_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, drntot'
     stop 1
  end if

  !**daily plant transpiration(npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'trans.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily rate of canopy transpiration (upper and lower) (mm/day)','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'adtrans','daily canopy transpiration','mm/day',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'adtrans',vardim,istart,icount,adtrans_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, adtrans'
     stop 1
  end if

  !**daily total evaporation (npoi,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'evap.nc'
  n3rd  = 1  !nlevel = 1
  if (new_file) then
    call inifile(idies,filen,'daily rate of evaporation (mm/day)','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'adevap','daily evaporation','mm/day',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'adevap',vardim,istart,icount,adevap_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, adevap'
     stop 1
  end if

  !**nitrate concentrations in solution(npoi,nsoilay,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'nconc.nc'
  n3rd      = nsoilay  !nlevel = nsoilay
  icount(3) = nsoilay
  vardim(2) = nsoilay
  long3rd   = 'soil depth'
  units3rd  = 'cm'
  if (new_file) then
    call inifile(idies,filen,'solute nitrate concentration','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'csoln','nitrate concentration each layer','layer',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'csoln',vardim,istart,icount,csoln_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, csoln'
     stop 1
  end if

  n3rd      = 1  !reset back to default
  icount(3) = 1
  vardim(2) = 1
  long3rd   = ''
  units3rd  = ''

  !**crop daily lai (npoi,npft,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'plai.nc'
  n3rd      = npft  !nlevel = npft
  icount(3) = npft
  vardim(2) = npft
  long3rd   = 'plant functional type'
  units3rd  = 'pft'
  if (new_file) then
    call inifile(idies,filen,'daily pft lai','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    !add global attribute to define pfts with text, use netcdf low-level command
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)

    call inivar(idies,'plai','plai for each pft','m2/m2',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'plai',vardim,istart,icount,plai_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, plai'
     stop 1
  end if

  n3rd      = 1  !reset back to default
  icount(3) = 1
  vardim(2) = 1
  long3rd   = ''
  units3rd  = ''
  
  !**total plant nitrogen uptake and stress factor by pft (npoi,npft,dyn_dim)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'plantn.nc'
  n3rd      = npft  !nlevel = npft
  icount(3) = npft
  vardim(2) = npft
  long3rd   = 'plant functional type'
  units3rd  = 'pft'
  if (new_file) then

    subfile = '-stressn'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    !add global attribute to define pfts with text, use netcdf low-level command
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'stressn'   ,'nitrogen stress factor for each pft'             ,'dimensionless',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-totnuptake'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'totnuptake','total nitrogen uptake for each pft'              ,'kg/m^2',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-fnleaf'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'fnleaf'    ,'nitrogen fraction in leaf for each  pft'         ,'fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-fnstem'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'fnstem'    ,'nitrogen fraction in stem for each  pft'         ,'fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-fnroot'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'fnroot'    ,'nitrogen fraction in root for each  pft'         ,'fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-fngrain'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'fngrain'   ,'nitrogen fraction in grain for each  pft'        ,'fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-fnplant'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'fnplant'   ,'nitrogen fraction in entire plant for each  pft' ,'fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)


    subfile = '-tnplant'
    call inifile(idies,trim(filen)//trim(subfile),'daily crop nitrogen variabiles','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)
    istat = NF_PUT_ATT_TEXT(idies,NF_GLOBAL,'pft_definition',npft*80,pftdef)
    call inivar(idies,'tnplant'   ,'nitrogen total in plant for each  pft'           ,'fraction',ndims,dimnames,OCEAN,istat)
    call endini(idies,istat)

  end if

  subfile = '-stressn'
  call write4dvar_ver2(trim(filen)//trim(subfile),'stressn'   ,vardim,istart,icount,stressn_dyn   ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-totnuptake'
  call write4dvar_ver2(trim(filen)//trim(subfile),'totnuptake',vardim,istart,icount,totnuptake_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-fnleaf'
  call write4dvar_ver2(trim(filen)//trim(subfile),'fnleaf'    ,vardim,istart,icount,fnleaf_dyn    ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-fnstem'
  call write4dvar_ver2(trim(filen)//trim(subfile),'fnstem'    ,vardim,istart,icount,fnstem_dyn    ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-fnroot'
  call write4dvar_ver2(trim(filen)//trim(subfile),'fnroot'    ,vardim,istart,icount,fnroot_dyn    ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-fngrain'
  call write4dvar_ver2(trim(filen)//trim(subfile),'fngrain'   ,vardim,istart,icount,fngrain_dyn   ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-fnplant'
  call write4dvar_ver2(trim(filen)//trim(subfile),'fnplant'   ,vardim,istart,icount,fnplant_dyn   ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  subfile = '-tnplant'
  call write4dvar_ver2(trim(filen)//trim(subfile),'tnplant'   ,vardim,istart,icount,tnplant_dyn   ,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, plantn.nc'
     stop 1
  end if

  n3rd      = 1  !reset back to default
  icount(3) = 1
  vardim(2) = 1
  long3rd   = ''
  units3rd  = ''
  
  !**upper and lower canopy daily lai (npoi,1,dyn_dim) for lai_dyn(:,1,:) and lai_dyn(:,2,:)
  filen = trim(myprocDir)//'output/daily/'//num//'_daily_'//'laicanopy.nc'
  n3rd      = 1  !nlevel = 1 for each variable
  icount(3) = 1
  vardim(2) = 1
  long3rd   = ''
  units3rd  = ''
  if (new_file) then
    call inifile(idies,filen,'daily lai-leaf area index of upper and lower vegetation canopies','ibis wdaily',cdate,nlonsub,lonscale,nlatsub,latscale, &
                 name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

    call inivar(idies,'laiu','daily lai of upper canopy','m2/m2',ndims,dimnames,OCEAN,istat)
    call inivar(idies,'lail','daily lai of lower canopy','m2/m2',ndims,dimnames,OCEAN,istat)

    call endini(idies,istat)
  end if

  call write4dvar_ver2(filen,'laiu',vardim,istart,icount,lai_dyn(:,1,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
  call write4dvar_ver2(filen,'lail',vardim,istart,icount,lai_dyn(:,2,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)

  if (istat .ne. 0) then
     write(*,*) 'ERROR in write_dailyout_dyn, laiu lail'
     stop 1
  end if

  n3rd      = 1  !reset back to default
  icount(3) = 1
  vardim(2) = 1
  long3rd   = ''
  units3rd  = ''

1131 continue

  !**daily4thmb variables - reset some variables for hourly4thmb
  layer_s = 5
  layer_e = 9
  if(flg_daily4thmb == 1 .or. flg_daily4thmb == 2) then

    if(flg_daily4thmb == 1) then  !daily-scale
      is_hourly = .false.
      n3rd      = 1

      filen = trim(myprocDir)//'output/daily/'//num//'_daily4thmb.nc'
    elseif(flg_daily4thmb == 2) then !hourly-scale
      is_hourly = .true.
      n3rd      = 24    !24 time steps per day
      vardim(2) = n3rd
      icount(3) = n3rd
      long3rd   = 'hourly'

      filen = trim(myprocDir)//'output/hourly/'//num//'_hourly4thmb.nc'
    end if

    if(new_file) then
      call inifile(idies,trim(filen),'Agro-IBIS daily input for THMB','ibis',cdate,nlonsub,lonscale,nlatsub,latscale, &
                   name3rd,long3rd,units3rd,n3rd,vals3rd(1:n3rd),pos3rd,tunits,'gregorian',istat)

      !**adsrunoff
      varName     = 'adsrunoff'
      longname    = 'daily average surface runoff'
      if(.not. is_hourly) then
        units     = 'mm/day'
      else
        units     = 'mm/hr'
      end if
      call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)

      !**addrainage
      varname     = 'addrainage'
      longname    = 'daily average drainage'
      if(.not. is_hourly) then
        units     = 'mm/day'
      else
        units     = 'mm/hr'
      end if
      call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)
      
      !do ilayer = 1, nsoilay
      do ilayer = layer_s, layer_e

        write(num_lay,'(i2)') ilayer
        varname     = 'addrainage_layer'//adjustl(num_lay)
        longname    = 'daily average drainage for layer '//adjustl(num_lay)
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)

      end do !on ilayer

      !**adrain
      varName     = 'adrain'
      longname    = 'daily average rain'
      if(.not. is_hourly) then
        units     = 'mm/day'
      else
        units     = 'mm/hr'
      end if
  
      call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)

      !**adevap
      varName     = 'adevap'
      longname    = 'daily average evaporation'
      if(.not. is_hourly) then
        units     = 'mm/day'
      else
        units     = 'mm/hr'
      end if
      call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)

      !**adtotnleach
      varName     = 'adtotnleach'
      longname    = 'daily average nitrogen leaching'
      if(.not. is_hourly) then
        units     = 'kg/ha/day'
      else
        units     = 'kg/ha/hr'
      end if
      call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)

      !do ilayer = 1, nsoilay
      do ilayer = layer_s, layer_e

        write(num_lay,'(i2)') ilayer
        varname     = 'adtotnleach_layer'//adjustl(num_lay)
        longname    = 'daily average nitrogen leaching for layer '//adjustl(num_lay)
        call inivar(idies,trim(varName),longname,units,ndims,dimnames,OCEAN,istat)

      end do !on ilayer

      call endini(idies,istat) !call endini after initializing ALL the variables (here only 1 variable is used)
    end if  !on if(new_file)

    !**write data (all vars here have dim (npoi,dyn_dim)

    !**adsrunoff
    varName = 'adsrunoff'

    if(is_hourly) then
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adsrunoff_hourly_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    else
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adsrunoff_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    end if
    if (istat .ne. 0) then
       write(*,*) 'ERROR in write_dailyout_dyn - daily4thmb, ',trim(varName)
       stop 1
    end if

    !**addrainage
    varName = 'addrainage'

    if(is_hourly) then
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,addrainage_layer_hourly_dyn(:,nsoilay,:,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    else
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,addrainage_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    end if
    if (istat .ne. 0) then
       write(*,*) 'ERROR in write_dailyout_dyn - daily4thmb, ',trim(varName)
       stop 1
    end if

    !do ilayer = 1, nsoilay
    do ilayer = layer_s, layer_e

      write(num_lay,'(i2)') ilayer
      varname = 'addrainage_layer'//adjustl(num_lay)

      if(is_hourly) then
        call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,addrainage_layer_hourly_dyn(:,ilayer,:,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
      else
        call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,addrainage_layer_dyn(:,ilayer,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
      end if

      if (istat .ne. 0) then
         write(*,*) 'ERROR in write_dailyout_dyn - daily4thmb, ',trim(varName)
         stop 1
      end if

    end do !on ilayer

    !**adrain
    varName = 'adrain'

    if(is_hourly) then
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adrain_hourly_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    else
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adrain_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    end if

    if (istat .ne. 0) then
       write(*,*) 'ERROR in write_dailyout_dyn - daily4thmb, ',trim(varName)
       stop 1
    end if

    !**adevap 
    varName = 'adevap'

    if(is_hourly) then
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adevap_hourly_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    else
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adevap_dyn,ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    end if

    if (istat .ne. 0) then
       write(*,*) 'ERROR in write_dailyout_dyn - daily4thmb, ',trim(varName)
       stop 1
    end if

    !**adtotnleach
    varName = 'adtotnleach'   !at default layer <isoilay> for compatability

    if(is_hourly) then
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adtotnleach_layer_hourly_dyn(:,isoilay,:,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    else
      call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adtotnleach_layer_dyn(:,isoilay,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
    end if

    if (istat .ne. 0) then
       write(*,*) 'ERROR in write_dailyout_dyn - daily4thmb, ',trim(varName)
       stop 1
    end if

    !do ilayer = 1, nsoilay
    do ilayer = layer_s, layer_e

      write(num_lay,'(i2)') ilayer
      varname = 'adtotnleach_layer'//adjustl(num_lay)

      if(is_hourly) then
        call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adtotnleach_layer_hourly_dyn(:,ilayer,:,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
      else
        call write4dvar_ver2(filen,trim(varName),vardim,istart,icount,adtotnleach_layer_dyn(:,ilayer,:),ftime(1:dyn_idx),tweight(1:dyn_idx),tdate,istat)
      end if

    end do !on ilayer

  end if !on if(flg_daily4thmb)
  !</write io completed>

  return
end subroutine write_dailyout_dyn
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : readnc_dailyout_freq
! purpose : read nc climate input for a month or year
!
! note: called only when flg_dailyout_freq != 0
!    
!
!                           Y.Li 2015-12-28
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine readnc_dailyout_freq(imonth,iyear,iwest,jnorth)

  use clim_file_utils
  use comgrid
  use compar
  use comwork,only: cdummy
  use cominput
  use combcs
  use daily_variables

  implicit NONE
  !--input
  integer:: imonth,iyear,iwest,jnorth
  !--local
  integer:: istart(4),icount(4),istart_time,icount_time,istat
  integer:: i,j,iyr,iyr_cal,ndays_pastMonths,inc

  character*200 filen        ! file name
  character*200 name3rd,aname

  integer(kind=4):: time_s,time_e          !clock counter, to check the total elapsed time
  integer(kind=4):: clock_rate, clock_max  !see system_clock() for reference

  logical first_time
  logical is_leap        ! Function: is the given year a leap year?

  save first_time

  call system_clock(time_s,clock_rate,clock_max)

  !--initialization
  first_time = .true.

  !year to be read info
  iyr = year_to_read(iyear, istyrd, overlap_start, overlap_end)

  if (iyr .lt. 0) then
    write(*,*) 'ERROR in readnc_dailyout_freq: tried to read a year prior to the start of the files'
    write(*,*) 'iyear = ', iyear
    write(*,*) 'istyrd = ', istyrd
    stop 1
  end if

  iyr_cal = iyr + istyrd
 
  !to ensure year to read cover days in iyear, i.e., when iyear is leap, so has iyr_cal
  inc = 1
  do while ( is_leap(iyear) .and. .not.is_leap(iyr_cal) )

    iyr_cal = iyr_cal + inc

    if(iyr_cal > overlap_end) then  !rewind back
      iyr_cal = istyrd
    end if

  end do

  if(first_time) then
    write(6,*) '(readnc_dailyout_freq)Reading daily climate data from year: ',iyr_cal
  end if !on if(first_time)

  if (iyear .lt. istyrd) then
    print *, 'daily data begins in year ', istyrd
    print *, 'not reading in daily data'
    return
  end if

  !get cummulative days before current month
  ndays_pastMonths = 0
  if(imonth > 1) then
    do i = 1, imonth-1
      ndays_pastMonths = ndays_pastMonths + ndaypm(i)
    end do
  end if

  if(flg_dailyout_freq == 1) then      !one month
    istart_time = ndays_pastMonths + 1
    icount_time = ndaypm(imonth)  !number of days for current month
  else if(flg_dailyout_freq == 2) then !one year
    istart_time = 1
    icount_time = ndaypy          !number of days for current year
  end if

  !dimension order (longitude,latitude,level,time)
  istart = (/iwest  ,jnorth ,1, istart_time /)
  icount = (/nlonsub,nlatsub,1, icount_time /)

  name3rd = trim(level_var)

  !**read daily precip
  if (fn_prec_daily .eq. ' ') then 
    write(*,*) 'ERROR in readnc_dailyout_freq: fn_prec_daily not specified'
    stop 1
  end if

  aname = var_prec_daily
  call climate_full_fname(filen, fn_prec_daily, iyr_cal, one_file_per_year)
  call readvar(filen,trim(aname),name3rd,istart,icount,xinprecd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
  if (istat .ne. 0) then
     write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_prec_daily'
     stop 1
  end if

  !**read tmin and tmax
  !set t and trng to un_init if this is the first time we're reading daily files (Note: unlike monthly variables, we do NOT calculate t and trng here, 
  !partly because of the difficulties that arise if these variables represent anomalies)
  if (read_tmin_tmax) then

    !read daily temp min
    if (fn_tmin_daily .eq. ' ') then
      write(*,*) 'ERROR in readnc_dailyout_freq: fn_tmin_daily must be specified for read_tmin_tmax = true'
      stop 1
    end if

    aname = var_tmin_daily
    call climate_full_fname(filen, fn_tmin_daily, iyr_cal, one_file_per_year)
    call readvar(filen,trim(aname),name3rd,istart,icount,xintmind_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
    if (istat .ne. 0) then
       write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_tmin_daily'
       stop 1
    end if

    !**read daily temp max
    if (fn_tmax_daily .eq. ' ') then
      write(*,*) 'ERROR in readnc_dailyout_freq: fn_tmax_daily must be specified for read_tmin_tmax = true'
      stop 1
    end if

    aname = var_tmax_daily
    call climate_full_fname(filen, fn_tmax_daily, iyr_cal, one_file_per_year)
    call readvar(filen,trim(aname),name3rd,istart,icount,xintmaxd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
    if (istat .ne. 0) then
       write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_tmax_daily'
       stop 1
    end if

    if (first_time) then
      call const(xintd,    npoi, un_init)
      call const(xintrngd, npoi, un_init)
    end if

  else ! .not. read_tmin_tmax
  !read t and trng; set tmin and tmax to un_init if this is the first time we're reading daily files (Note: unlike monthly variables, we do NOT calculate tmin and tmax here, 
  !partly because of the difficulties that arise if these variables represent anomalies)

    !** read daily temp 
    if (fn_temp_daily .eq. ' ') then
      write(*,*) 'ERROR in readnc_dailyout_freq: fn_temp_daily must be specified for read_tmin_tmax = false'
      stop 1
    end if

    aname = var_temp_daily
    call climate_full_fname(filen, fn_temp_daily, iyr_cal, one_file_per_year)
    call readvar(filen,trim(aname),name3rd,istart,icount,xintd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
    if (istat .ne. 0) then
       write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_temp_daily'
       stop 1
    end if

    !**read daily trange
    if (fn_dtr_daily .eq. ' ') then
      write(*,*) 'ERROR in readnc_dailyout_freq: fn_dtr_daily must be specified for read_tmin_tmax = false'
      stop 1
    end if

    aname = var_dtr_daily
    call climate_full_fname(filen, fn_dtr_daily, iyr_cal, one_file_per_year)
    call readvar(filen,trim(aname),name3rd,istart,icount,xintrngd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
    if (istat .ne. 0) then
       write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_dtr_daily'
       stop 1
    end if

    if (first_time) then
      call const(xintmind, npoi, un_init)
      call const(xintmaxd, npoi, un_init)
    end if

  end if !on if read_tmin_tmax

  !read rads, set cloud to un_init if this is the first time we're reading daily files
  if (read_radiation) then

    !**read daily radiation
    if (fn_rads_daily .eq. ' ') then
      write(*,*) 'ERROR in readnc_dailyout_freq: fn_rads_daily must be specified for read_radiation = true'
      stop 1
    end if

    aname = var_rads_daily
    call climate_full_fname(filen, fn_rads_daily, iyr_cal, one_file_per_year)
    call readvar(filen,trim(aname),name3rd,istart,icount,xinradsd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
    if (istat .ne. 0) then
       write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_rads_daily'
       stop 1
    end if

    if (first_time) then
      call const(xincldd, npoi, un_init)
    end if

  else ! .not. read_radiation
  !Read cloud, set rads to un_init if this is the first time we're reading daily files

    !**read daily cloudiness
    if (fn_cloud_daily .eq. ' ') then
      write(*,*) 'ERROR in readnc_dailyout_freq: fn_cloud_daily must be specified for read_radiation = false'
      stop 1
    end if

    aname = var_cloud_daily
    call climate_full_fname(filen, fn_cloud_daily, iyr_cal, one_file_per_year)
    call readvar(filen,trim(aname),name3rd,istart,icount,xincldd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
    if (istat .ne. 0) then
       write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_cloud_daily'
       stop 1
    end if

    if (first_time) then
      call const(xinradsd, npoi, un_init)
    end if

  end if !on if read_radiation

  !**read daily windspeed
  if (fn_wspd_daily .eq. ' ') then
    write(*,*) 'ERROR in readnc_dailyout_freq: fn_wspd_daily not specified'
    stop 1
  end if

  aname = var_wspd_daily
  call climate_full_fname(filen, fn_wspd_daily, iyr_cal, one_file_per_year)
  call readvar(filen,trim(aname),name3rd,istart,icount,xinwindd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
  if (istat .ne. 0) then
     write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_wspd_daily'
     stop 1
  end if

  !**read daily humidity
  if (fn_rh_daily .eq. ' ') then
    write(*,*) 'ERROR in readnc_dailyout_freq: fn_rh_daily not specified'
    stop 1
  end if

  aname = var_rh_daily
  call climate_full_fname(filen, fn_rh_daily, iyr_cal, one_file_per_year)
  call readvar(filen,trim(aname),name3rd,istart,icount,xinqd_dyn,cdummy(1),cdummy(nlonsub+1),cdummy(2*nlonsub+1),cdummy(3*nlonsub+1),istat)
  if (istat .ne. 0) then
     write(*,*) '(readnc_dailyout_freq)ERROR in readvar, fn_rh_daily'
     stop 1
  end if


  !all read-in completed
  first_time = .false.
  call system_clock(time_e,clock_rate,clock_max)
  write(6,'(A,g16.8,A)') '(readnc_dailyout_freq)time to read chunk daily climate data ',real(time_e-time_s)/(real(clock_rate)),' seconds'

  return
end subroutine readnc_dailyout_freq
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : assign_dailyout_freq
! purpose : assign read-in chunk data (a month or year) to daily basis
!
! note: called only when flg_dailyout_freq != 0
!    
!
!                           Y.Li 2015-12-29
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine assign_dailyin_freq

  use cominput
  use comgrid
  use combcs
  use compar
  use daily_variables

  implicit NONE

  !--begin the work

  call arr2vec( xinprecd_dyn(:,dyn_idx)  , xinprecd(1) )   !**daily precip

  if (read_tmin_tmax) then
    call arr2vec( xintmind_dyn(:,dyn_idx), xintmind(1) )   !**daily temp min
    call arr2vec( xintmaxd_dyn(:,dyn_idx), xintmaxd(1) )   !**daily temp max    
  else !.not. read_tmin_tmax
    call arr2vec( xintd_dyn(:,dyn_idx)   , xintd(1) )      !**daily temp
    call arr2vec( xintrngd_dyn(:,dyn_idx), xintrngd(1) )   !**daily trange
  end if !on if read_tmin_tmax

  if (read_radiation) then
    call arr2vec( xinradsd_dyn(:,dyn_idx), xinradsd(1) )   !**daily radiation
  else !.not. read_radiation
    call arr2vec( xincldd_dyn(:,dyn_idx), xincldd(1) )     !**daily cloudiness
  end if !on if read_radiation

  call arr2vec( xinwindd_dyn(:,dyn_idx), xinwindd(1) )     !**daily windspeed
  call arr2vec( xinqd_dyn(:,dyn_idx), xinqd(1) )           !**daily humidity

  return
end subroutine assign_dailyin_freq
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : preproc_management
! purpose : pre-processing for management
!
! note: 
!
!                           Y.Li 2016-04-11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine preproc_management
  use compar,only: flg_mngt

  call init_management_variables

  if(flg_mngt == 1) then
    write(6,'(/,1x,A)') 'Using management control flg_mngt == 1'
    call read_management_nml()
    call alloc_management_variables()
    call processing_management()
  end if

  return
end subroutine preproc_management
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : init_management_variables
! purpose : initialize management variables
!
! note    : 
!                           Y.Li 2016-04-11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine init_management_variables
  use management_para

  implicit NONE
  
  !namelist-climate
  climate_mngt_mode    = 'none'
  climate_mngt_file    = 'none'

  !namelist-climate
  fertilizer_mngt_mode       = 'none'
  fertilizer_mngt_year_start = 1945
  fertilizer_mngt_year_end   = 1945
  fertilizer_mngt_file       = 'none'

  !climate variables
  climateMngt_yearsList = 9999

  !fertilizer variables
  fertMngt_ys_idx = 0
  fertMngt_ye_idx = 0
  fertMngt_nyears = 0

  return
end subroutine init_management_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : read_management_nml
! purpose : read namelist <management.nml>
!
! note    : 
!                           Y.Li 2016-04-11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine read_management_nml
  use management_para

  implicit NONE
  !--local
  integer::     iunit,ierr
  character*200 file_nml,argc

  !--define namelists
  namelist/climate/climate_mngt_mode,climate_mngt_file
  namelist/fertilizer/fertilizer_mngt_mode,fertilizer_mngt_year_start,fertilizer_mngt_year_end,fertilizer_mngt_file

  !--read namelist <management.nml>
  file_nml = 'management.nml'  !default value

  write(6,*) 'reading management namelist file <management.nml>'

  iunit=8
  open(unit=iunit,file=trim(file_nml),status='old',action='read',iostat=ierr)
  if(ierr /= 0) then
     write(6,*) '(read_management_nml)Error opening ',trim(file_nml)
     write(6,*) 'No such file'
     stop 'Stopping the program...'
  endif

  read(iunit,nml=climate)
  rewind(iunit)

  read(iunit,nml=fertilizer)
  rewind(iunit)

  close(iunit)

  return
end subroutine read_management_nml
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : processing_management
! purpose : process management based on read-in <management.nml>
!
! note    : called only when flg_mngt = 1
!
!                           Y.Li 2016-04-11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine processing_management

  implicit NONE

  call climate_management()
  call fertilizer_management()


end subroutine processing_management
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : climate_management
! purpose : 
!
! note    : 
! 1) climate_mngt_mode = 'dat' format:
!    <nyears> 
!    <iyear>    <climate_year>
!    <iyear>    <climate_year>
! prerequisite: (the code will check and stop if not)
!     - <nyears>      : greater or equal to <nrun> 
!     - <iyear>       : follow the same simulation years, in order
!     - <climate_year>: if<iyear> is leap year, <climate_year> has to be leap year 
!
!                           Y.Li 2016-04-11
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine climate_management
  use management_para
  use compar,only: irestart_glo,irstyear_glo,iyear0_glo,nrun_glo

  implicit NONE
  !--local
  integer:: i,tmp_iyear,sim_year,iunit,ierr,mngt_nyears
  logical:: is_leap

  iunit = 8 
  if(trim(climate_mngt_mode) == 'dat') then
    open(unit=iunit,file=trim(climate_mngt_file),status='old',action='read',iostat=ierr)
    if(ierr /= 0) then
       write(6,*) '(climate_management)Error opening ',trim(climate_mngt_file)
       write(6,*) 'No such file'
       stop 'Stopping the program...'
    endif

    if(irestart_glo == 1) then
      sim_year = irstyear_glo
    else
      sim_year = iyear0_glo
    end if

    read(iunit,*) mngt_nyears
    if(mngt_nyears < nrun_glo) then
      write(6,'(A)') '(climate_management)Error! 1st line: <nyears> less than <nrun>, need more management years'
      write(6,*) 'Stopping the program...'
      stop
    end if 

    if(mngt_nyears > MAX_YEARS) then
      write(6,'(A)') '(climate_management)Error! number of years exceed MAX_YEARS defined. Increase MAX_YEARS first'
      write(6,'(A,i4,A,i4)') 'number of years: ',mngt_nyears,'MAX_YEARS: ',MAX_YEARS 
      write(6,*) 'Stopping the program...'
      stop
    end if
      
    do i = 1, nrun_glo   !use nrun_glo as we just need that many (mngt_nyears>=nrun)
      read(iunit,*) tmp_iyear,climateMngt_yearsList(i)

      if(tmp_iyear /= sim_year) then
        write(6,*) '(climate_management)Error! at ith year: ',i
        write(6,'(A,i4,A,i4)') 'Read-in year = ',tmp_iyear, ' not consistent with simulation year = ',sim_year
        write(6,*) 'Stopping the program...'
        stop
      end if

      if(is_leap(tmp_iyear)) then
        if(.not. is_leap(climateMngt_yearsList(i))) then
          write(6,*) '(climate_management)Error! at ith year ',i
          write(6,'(A,i4,A,i4,A)') 'simulation year = ',tmp_iyear, ' is leap year, climate year = ',climateMngt_yearsList(i), ' is not...'
          write(6,*) 'Stopping the program...'
          stop
        end if
      end if

      sim_year = sim_year+1
    end do  !on do i =1, mngt_nyears
    
    close(iunit)

    write(6,*) 'parsing climate management file: ',trim(climate_mngt_file),' completed...'
  end if !on if trim(climate_mngt_mode) == 'dat'

  
  return
end subroutine climate_management
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : fertilizer_management
! purpose : 
!
! note    : 
!
!                           Y.Li 2016-04-20
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine fertilizer_management
  use management_para

  implicit NONE
  !--local
  integer:: iunit,ierr,tmp_year_start,tmp_nyears
  integer:: i,nlines,iline(367)
  real   :: tmp_sum
  character*200 filen, tmp_buf

  iunit = 8 

  if(trim(fertilizer_mngt_mode) == 'years_ave') then

    !collect basic info
    filen = 'input/frate.hdr'
    open(unit=iunit,file=trim(filen),status='old',action='read',iostat=ierr)
    if(ierr /= 0) then
       write(6,*) '(fertilizer_management)Error opening ',trim(filen)
       write(6,*) 'No such file'
       stop 'Stopping the program...'
    endif

    read(iunit,*) tmp_year_start
    read(iunit,*) tmp_nyears

    close(iunit)

    !checking
    if(fertilizer_mngt_year_start < tmp_year_start .or. fertilizer_mngt_year_start > (tmp_year_start+tmp_nyears+1) )then 
      write(6,'(A,i4,A)') 'fertilizer_management)Error! <fertilizer_mngt_year_start> = ',fertilizer_mngt_year_start, ', not within available range'
      write(6,*) 'Stopping the program...'
      stop
    end if

    if(fertilizer_mngt_year_end < tmp_year_start .or. fertilizer_mngt_year_end > (tmp_year_start+tmp_nyears+1) )then 
      write(6,'(A,i4,A)') 'fertilizer_management)Error! <fertilizer_mngt_year_end> = ',fertilizer_mngt_year_end, ', not within available range'
      write(6,*) 'Stopping the program...'
      stop
    end if

    if(fertilizer_mngt_year_end < fertilizer_mngt_year_start) then
      write(6,'(A)') 'fertilizer_management)Error! <fertilizer_mngt_year_end> less than <fertilizer_mngt_year_start>'
      write(6,*) 'Stopping the program...'
      stop
    end if

    fertMngt_ys_idx = fertilizer_mngt_year_start - tmp_year_start + 1
    fertMngt_ye_idx = fertilizer_mngt_year_end   - tmp_year_start + 1
    fertMngt_nyears = fertMngt_ye_idx-fertMngt_ys_idx+1

    write(6,'(A,i4,2x,i4)') 'Using fertilizer management: average value between years ',fertilizer_mngt_year_start, fertilizer_mngt_year_end

  end if !on if fertilizer_mngt_mode == 'years_ave'

  if(trim(fertilizer_mngt_mode) == 'dat') then

    write(6,'(A,A)') 'Using fertilizer management: daily distribution since planting date by file: ',trim(fertilizer_mngt_file)

    open(unit=iunit,file=trim(fertilizer_mngt_file),status='old',action='read',iostat=ierr)
    if(ierr /= 0) then
       write(6,*) '(fertilizer_management)Error opening ',trim(fertilizer_mngt_file)
       write(6,*) 'No such file'
       write(6,*) 'Stopping the program...'
       stop
    endif

    !read the header; following 367 lines of data
    read(iunit,*) 

    tmp_sum = 0.0
    nlines = 367
    do i = 1, nlines
      read(iunit,*,iostat=ierr) iline(i),fertMngt_pct(i)

      if(ierr /= 0) then
       write(6,*) '(fertilizer_management)Error reading file ',trim(fertilizer_mngt_file)
       write(6,*) 'at line# ',i+1
       write(6,*) 'Stopping the program...'
       stop
      end if 

      !checking
      if( iline(i)+1-i /=0 ) then
        write(6,*) '(fertilizer_management)Error at line # ',i+1
        write(6,'(1x,A,i4,A,i4)') 'days since planting date not consistent. current value = ',iline(i), ', expecting = ',i-1
        write(6,*) 'Stopping the program...'
        stop
      end if

      tmp_sum = tmp_sum + fertMngt_pct(i)

    end do !on i = 1, nlines

    close(iunit)

    !checking
    if(abs(tmp_sum-1.0) > 1e-10) then
       write(6,*) '(fertilizer_management)Error'
       write(6,'(A,g12.5)') 'Summation of <fertMngt_pct> is not equal to 1. Current sum = ',tmp_sum
       stop 'Stopping the program...'
    end if 


    write(6,*) 'parsing fertilizer management file: ',trim(fertilizer_mngt_file),' completed...'

  end if !on if fertilizer_mngt_mode == 'dat'

  return
end subroutine fertilizer_management
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : alloc_management_variables
! purpose : 
!
! note    : 
!
!                           Y.Li 2016-04-20
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine alloc_management_variables
  use management_para
  use comgrid,only: npoi
  use compar,only: flg_mngt

  implicit NONE

  !--local
  integer:: ierr

  !fertilizer variables
  if(flg_mngt == 1 .and. trim(fertilizer_mngt_mode) == 'dat') then
    allocate(mngt_fertnitro0_13(npoi),STAT=ierr)
    allocate(mngt_fertnitro0_14(npoi),STAT=ierr)
    allocate(mngt_fertnitro0_15(npoi),STAT=ierr)

    allocate(fertMngt_pct(367),STAT=ierr)
  else
    allocate(mngt_fertnitro0_13(1),STAT=ierr)
    allocate(mngt_fertnitro0_14(1),STAT=ierr)
    allocate(mngt_fertnitro0_15(1),STAT=ierr)

    allocate(fertMngt_pct(1),STAT=ierr)
  end if !on if flg_mngt & fertilizer_mngt_mode


  return
end subroutine alloc_management_variables
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! name    : dealloc_management_variables
! purpose : 
!
! note    : 
!
!                           Y.Li 2016-04-20
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine dealloc_management_variables
  use management_para

  implicit NONE
  !--local
  integer:: ierr

  !fertilizer variables
  deallocate(mngt_fertnitro0_13, STAT=ierr)
  deallocate(mngt_fertnitro0_14, STAT=ierr)
  deallocate(mngt_fertnitro0_15, STAT=ierr)

  return
end subroutine dealloc_management_variables
