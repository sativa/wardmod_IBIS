c Bill Sacks, Jan. 2010

c----------------------------------------------------------
      subroutine read_input_description (datasource)
c----------------------------------------------------------
c
c Read a namelist describing the input data source.
c The namelist is contained in the file 'input_descriptions/<datasource>'
c
c uses:
c
      use clim_file_utils
      use cominput
      use compar,only: nc_file_bufrsize_KB
c
      implicit none
c
c Arguments
c
      character*256 datasource
c
c Local variables
c
      logical verbose ! if true, print out extra info
      integer funit   ! file unit assignment for input
      character*64 input_directory  ! directory containing input descriptions
      integer ios     ! I/O error status

      parameter (verbose = .FALSE.)
      parameter (funit = 101)
      parameter (input_directory = 'input_descriptions')
c 
c Namelist for input
c 
      namelist /input_description/ read_daily_directly,
     >     one_file_per_year, read_tmin_tmax, read_radiation,
     >     correct_sphumid, conserve_tdew, istyrm, nanom, istyrd,
     >     nanomd, level_var,
     >     
     >     fn_cloud_mon_clim, fn_cloud_mon, fn_cloud_daily,
     >     fn_rads_mon_clim, fn_rads_mon, fn_rads_daily,
     >     fn_temp_mon_clim, fn_temp_mon, fn_temp_daily,fn_dtr_mon_clim,
     >     fn_dtr_mon, fn_dtr_daily, fn_tmin_mon_clim, fn_tmin_mon,
     >     fn_tmin_daily, fn_tmax_mon_clim, fn_tmax_mon,fn_tmax_daily,
     >     fn_prec_mon_clim, fn_prec_mon, fn_prec_daily, fn_rh_mon_clim,
     >     fn_rh_mon, fn_rh_daily, fn_wspd_mon_clim,fn_wspd_mon,
     >     fn_wspd_daily, fn_wetd_mon_clim, fn_wetd_mon, fn_tminavgann,
     >     
     >     var_cloud_mon_clim, var_cloud_mon, var_cloud_daily,
     >     var_rads_mon_clim, var_rads_mon, var_rads_daily,
     >     var_temp_mon_clim, var_temp_mon,
     >     var_temp_daily,var_dtr_mon_clim, var_dtr_mon, var_dtr_daily,
     >     var_tmin_mon_clim, var_tmin_mon, var_tmin_daily,
     >     var_tmax_mon_clim, var_tmax_mon,var_tmax_daily,
     >     var_prec_mon_clim, var_prec_mon, var_prec_daily,
     >     var_rh_mon_clim, var_rh_mon, var_rh_daily,
     >     var_wspd_mon_clim,var_wspd_mon, var_wspd_daily,
     >     var_wetd_mon_clim, var_wetd_mon, var_tminavgann,
     >     nc_file_bufrsize_KB

c ---------------------------------------------------------------
c Set default values
c ---------------------------------------------------------------

      read_daily_directly = .FALSE.
      one_file_per_year = .FALSE.
      read_tmin_tmax = .FALSE.
      read_radiation = .FALSE.
      correct_sphumid = .FALSE.
      conserve_tdew = .FALSE.

      istyrm = 9999
      nanom = 0
      istyrd = 9999
      nanomd = 0

      level_var = 'level'

      fn_cloud_mon_clim = ' '
      fn_cloud_mon = ' '
      fn_cloud_daily = ' '
      fn_rads_mon_clim = ' '
      fn_rads_mon = ' '
      fn_rads_daily = ' '
      fn_temp_mon_clim = ' '
      fn_temp_mon = ' '
      fn_temp_daily = ' '
      fn_dtr_mon_clim = ' '
      fn_dtr_mon = ' '
      fn_dtr_daily = ' '
      fn_tmin_mon_clim = ' '
      fn_tmin_mon = ' '
      fn_tmin_daily = ' '
      fn_tmax_mon_clim = ' '
      fn_tmax_mon = ' '
      fn_tmax_daily = ' '
      fn_prec_mon_clim = ' '
      fn_prec_mon = ' '
      fn_prec_daily = ' '
      fn_rh_mon_clim = ' '
      fn_rh_mon = ' '
      fn_rh_daily = ' '
      fn_wspd_mon_clim = ' '
      fn_wspd_mon = ' '
      fn_wspd_daily = ' '
      fn_wetd_mon_clim = ' '
      fn_wetd_mon = ' '
      fn_tminavgann = ' '

      var_cloud_mon_clim = 'cld'
      var_cloud_mon = 'cld'
      var_cloud_daily = 'cld'
      var_rads_mon_clim = 'rads'
      var_rads_mon = 'rads'
      var_rads_daily = 'rads'
      var_temp_mon_clim = 'temp'
      var_temp_mon = 'temp'
      var_temp_daily = 'temp'
      var_dtr_mon_clim = 'dtr'
      var_dtr_mon = 'dtr'
      var_dtr_daily = 'dtr'
      var_tmin_mon_clim = 'tmin'
      var_tmin_mon = 'tmin'
      var_tmin_daily = 'tmin'
      var_tmax_mon_clim = 'tmax'
      var_tmax_mon = 'tmax'
      var_tmax_daily = 'tmax'
      var_prec_mon_clim = 'prec'
      var_prec_mon = 'prec'
      var_prec_daily = 'prec'
      var_rh_mon_clim = 'rh'
      var_rh_mon = 'rh'
      var_rh_daily = 'rh'
      var_wspd_mon_clim = 'wspd'
      var_wspd_mon = 'wspd'
      var_wspd_daily = 'wspd'
      var_wetd_mon_clim = 'wetd'
      var_wetd_mon = 'wetd'
      var_tminavgann = 'tminavgann'

      nc_file_bufrsize_KB = 8  !in KB 

c ---------------------------------------------------------------
c Open file and read it
c ---------------------------------------------------------------

      open(unit=funit, 
     >     file=trim(input_directory) // '/' // trim(datasource), 
     >     status='old',
     >     iostat=ios)

      if (ios .ne. 0) then
        write(*,*) 'READ_INPUT_DESCRIPTION: Error opening file: ',
     >       trim(input_directory) // '/' // trim(datasource)
        write(*,*) 'IOS = ', ios
        stop 1
      end if

      read(funit, nml=input_description, iostat=ios)

      if (ios .ne. 0) then
        write(*,*) 'READ_INPUT_DESCRIPTION: Error reading file: ',
     >       trim(input_directory) // '/' // trim(datasource)
        write(*,*) 'IOS = ', ios
        stop 1
      end if

      close(funit)

c ---------------------------------------------------------------
c Do some consistency checks
c ---------------------------------------------------------------

      if (read_daily_directly .and. nanom .gt. 0) then
        write(*,*) 'ERROR DETECTED IN READ_INPUT_DESCRIPTION:'
        write(*,*) 'For read_daily_directly = true, nanom should be 0'
        stop 1
      end if

c ---------------------------------------------------------------
c Set auxiliary variables
c ---------------------------------------------------------------

      call overlap(overlap_start, overlap_end, istyrm, nanom,
     >     istyrd, nanomd)

      call climate_years_avail(avail_leaps, avail_nonleaps,
     >     overlap_start, overlap_end)
      
c ---------------------------------------------------------------
c Print values, if desired
c ---------------------------------------------------------------

      if (verbose) then
        write(*, nml=input_description)
        write(*,*)
        write(*,*) 'overlap_start = ', overlap_start
        write(*,*) 'overlap_end = ', overlap_end
        write(*,*) 'avail_leaps = '
        write(*,*) avail_leaps
        write(*,*) 'avail_nonleaps = '
        write(*,*) avail_nonleaps
      end if

      return
      end
