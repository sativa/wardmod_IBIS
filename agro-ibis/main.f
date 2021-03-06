1844167voc
c #    #    ##       #    #    #
c ##  ##   #  #      #    ##   #
c # ## #  #    #     #    # #  #
c #    #  ######     #    #  # #
c #    #  #    #     #    #   ##
c #    #  #    #     #    #    #
c
c ---------------------------------------------------------------
      program main
c ---------------------------------------------------------------
c
c uses:
c 
      use comgrid
      use compar
      use cominput
      use comatm
      use comdiag
      use comveg
      use comcrop
      use comnitr
      use comsum,only:adsrunoff,addrainage,adrain,adevap
      use combcs,only:cmask
      use management_para

! S. Slavin added following module in May 2016 for read optimization.
! flg_rdopt is declared in compar.f and is also used in ies-io.f
      use rdopt

      implicit none
c
c local variables
c
      integer icount,        ! number of times ran2 is called (for restart)
     >        icountdum,     !
     >        istep,         ! timestep counter (per day)
     >        iday,          ! daily loop counter
     >        imonth,        ! monthly loop counter
     >        iwest,         ! 1st lon index for subset
     >        iy1,           ! first year for year loop
     >        iy2,           ! last year for year loop
     >        iyear,         ! yearly loop counter
     >        iyear_climate, ! for spinup with weather_generator=FALSE: year of climate data to use
     >        iyear0,        ! first year of simulation
     >        iyranom,       ! year to start reading monthly anomalies (ignored if weather_generator=FALSE)
     >        iyrdaily,      ! year to start reading daily means
     >        iyrlast,       ! last year of previous run (for restart)
     >        idiag,         ! number of diagnostic files requested
     >        idailyout,     ! 0: no daily output 1: daily output
     >        imonthout,     ! 0: no monthly output 1: monthly output
     >        flg_wrestart,  ! 1: write restart files (default); 0: do not write restart files
     >        iyearout,      ! 0: no yearly output 1: yearly output
     >        irestart,      ! 0: normal mode 1: restart mode
     >        irstyear,      ! actual calendar year that restart occurs
     >        iholdsoiln,    ! 0: don't save inorganic soil N values 1: keep inorganic soil N values 
     >        isimveg,       ! 0: static veg 1: dynam veg initialized w/ fixed
     >                       ! 2: dynam veg initialized w/ cold start
     >        isimfire,      ! 0: fixed fire  1: dynamic fire
     >        isimco2,       ! 0: fixed co2   1: changing co2
     >        jday,          ! julian day of the simulation
     >        jnorth,        ! 1st lat index for subset
     >        niter,         ! total number of time iterations per day
     >        nday,          ! number of days since the start of the simulation
     >        ndayr,         !
     >        nrun,          ! # of years in this run
     >        nspinsoil,     ! year of simulation that soil c reaches equilibrium 
     >        plen,          ! length of precipitation event in timesteps (see plens)
     >        plenmax,       ! upper bound on plen
     >        plenmin,       ! lower bound on plen
     >        seed,          ! value used to initialize the random # generator
     >        spin,          ! counter for iterations used to spin up soilbgc
     >        spinmax,       ! maximum iterations used to spin up soilbgc
     >        eqyears,       !
     >        soilcspin,     ! 0: no spinup procedure for soil c  1: acceleration procedure used
     >        lun, ifile,    ! file indices
     >        i, n, ivar,    ! loop indices
     >        cropsums,
     >        nerr           ! number of errors found in user input
c
      real    co2init,       ! initial co2 concentration in mol/mol
     >        o2init,        ! initial o2 concentration in mol/mol
     >        dran,          ! random number from generator
     >        plens,         ! length of precipitation event in seconds
     >        startp,        ! time to start precipitation event (seconds since midnight)
     >        endp,          ! time to end precipitation event (seconds since midnight)
     >        ilens,         ! length of irrigation event in seconds
     >        starti,        ! time to start irrigation event (seconds since midnight)
     >        endi,          ! time to end irrigation event (seconds since midnight) 
     >        spincons,      ! # times soil gbc called each day during spinup
     >        spinfrac,      ! fraction of nspinsoil time with max iteration spinmax
     >        slope,         ! rate of decrease of number of iteration for soil C spinup 
     >        snorth,        ! north latitude for subsetting std grid
     >        ssouth,        ! south latitude for subsetting std grid
     >        swest,         ! west longitude for subsetting std grid
     >        seast,         ! east longitude for subsetting std grid
     >        ffact,         ! numerical multiplying factor applied to N fertilizer being applied after end of fertilizer data
     >        test,          ! test on timestep 
     >        time           ! time of day since midnight (in seconds)

c should we use the weather generator for spin-up? (if false, then in
c each year prior to iyrdaily, we choose a random year of available
c daily climate data [but only select from leap years if this is a leap
c year, and vice versa], and use that year's data)
      logical weather_generator
c
c External
c
      integer random_seed    ! Function: get a seed for the random number generator based on CPU time
      
      real ran2,             ! Function : random number generator
     >     get_time          ! Function: get current time

      logical is_leap        ! Function: is the given year a leap year?
c
      integer imetyear,      ! year to start using hourly met data for climate 
     >        imetend        ! ending year of reading in hourly met data
c
c file in input_descriptions directory describing the climate data source
      character*256 datasource

      character*32 
     >  planting_input,  ! Suffix on the planting date input file name (e.g., 'transient', 'detrended', 'linearized' or 'fixed')
     >  cultivar_input   ! Suffix on the cultivar input file name (e.g., 'transient' or 'fixed')

      character*39 pathn
      character*80 filen

c variables for performance checking (Y.Li)
      character*8         :: p_date
      character*10        :: P_time
      character*5         :: p_zone
      integer,dimension(8):: p_values

      integer(kind=8)     :: clock_s,clock_e        !clock counter, to check elapsed time for section of interest
      integer(kind=8)     :: clock_ts,clock_te      !clock counter, to check the total elapsed time for the whole program
      integer(kind=8)     :: yearly_ts,yearly_te    !clock counter, to check the total elapsed time for a year
      integer(kind=8)     :: clock_rate, clock_max  !see system_clock() for reference

c other local variables (Y.Li)
      integer:: year_idx   !from 1 to nrun
      integer:: length

      integer:: nargc,iarg
      character*200 fnml,argc   !file name of namelist (default='ibis.infile')
      character*200 str_buf     !temporary string buffer

      logical:: flg_read_daily  !read daily climate input in a monthly or yearly frequency

      !for flg_ghost_io
      character*200  myprocDir_ori  !orginal path before prefixing the ghost_io
      character*200 message         !for system command
      integer:: exit_s, cmd_s 
      character*4 num               
c flg_sys_stats Added by Cicada Dennis, May 2016.
      integer:: flg_sys_stats   !Report disk usage each year? 1=yes, 0=no (default)      

c number of days per month
c
      data snorth, ssouth, swest, seast / 90., -90., -180., 180. /
      ndaypm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      pathn = '/Volumes/wet/kucharik/IBIS_input/input/'

c Initialize inputs
      flg_wrestart = 1

      isoil_layer  = 0
      iyr_pft      = 0
      iyr_fert     = 0
      ihourlyout   = 0
    
      iagro_thmb     = 0
      agro_thmb_path = './'
      commu_MaxWait  = 200   
      commu_sleep    = 1

      checked_grid   = .false.
      error_showed_sumday = .false.
      error_showed_plen   = .false.

      flg_udg = 0
      udg_is  = 0
      udg_ie  = 0
      udg_js  = 0
      udg_je  = 0

      lon_llcorner = 0.0
      lat_llcorner = 0.0

      flg_debug         = 0
      debug_ncounts     = 1
      flg_detail        = 0
      flg_dailyin_freq  = 0
      flg_dailyout_freq = 0
      flg_daily4thmb    = 0
      flg_sandlayers    = 0
      flg_cmask         = 0
      flg_coverCrops    = 0
      flg_mngt          = 0
 
      flg_mpi           = 0  !single processor run

      flg_read_daily    = .false.

      call system_clock(clock_ts,clock_rate,clock_max)

c
c ---------------------------------------------------------------
c                       j o b   c o n t r o l 
c ---------------------------------------------------------------
c
c open local input file called 'ibis.infile'
c and read in pertinent information
c
      !Enable readin commandline argument for input file <ibis.infile>
      !for parallelization (Y.Li) & implicit IO operations for Big Red 2 (ghost_io)

      fnml='ibis.infile'    !default
      flg_ghost_io = 0      !default, normal IO
      ! Flag added by Shawn Slavin.
      flg_rdopt         = 0 !default, does not use read optimization
      ! Flag added by Cicada Dennis.
      flg_sys_stats     = 0 !default, does not report system statistics.

      nargc = command_argument_count()
      do iarg = 1, nargc

        call get_command_argument(iarg,argc)
        if(trim(argc)=='-f') then
          call get_command_argument(iarg+1,argc)
          read(argc,'(A)') fnml
          write(6,*) 'read-in infile: ',trim(fnml)
        end if  !on if(trim(argc)

        if(trim(argc) =='--ghost_io') then
          flg_ghost_io = 1
          fnml = '/tmp/'//fnml
        end if

        ! S. Slavin added --rdopt flag.
        if(trim(argc) == '--rdopt') then
           flg_rdopt = 1
        end if

        ! Cicada Dennis added --sys_stats flag.
        if(trim(argc) =='--sys_stats') then
          flg_sys_stats = 1
        end if

      end do
     
      ! S. Slavin added if statement.
      if ( flg_rdopt == 0 ) then
         write(*,'(a)') 'Read optimization for readvar() set OFF'
      else
         write(*,'(a)') 'Read optimization for readvar() set ON'
      end if
     
      lun = 12
      write(0,*) 'About to open ibis.infile for reading: ',trim(fnml)
      write(6,*) 'About to open ibis.infile for reading: ',trim(fnml)
      open (lun, status='old', file=trim(fnml))
c
      read (lun,*) irestart
      read (lun,*) irstyear
      read (lun,*) iyear0
      read (lun,*) nrun
      read (lun,*) datasource
      read (lun,*) weather_generator
      read (lun,*) iyranom
      read (lun,*) iyrdaily
      read (lun,*) imetyear
      read (lun,*) imetend
      read (lun,*) seed
      read (lun,*) soilcspin
      read (lun,*) flg_wrestart
      read (lun,*) iyearout
      read (lun,*) imonthout
      read (lun,*) idailyout
      read (lun,*) ihourlyout
      read (lun,*) isimveg
      read (lun,*) isimfire
      read (lun,*) isimco2
      read (lun,*) irrigate_flag
      read (lun,*) nstress
      read (lun,*) isoybean
      read (lun,*) imaize
      read (lun,*) iwheat
      read (lun,*) irotation
      read (lun,*) iholdsoiln
      read (lun,*) isoil_layer   
      read (lun,*) iyr_pft       
      read (lun,*) iyr_fert      
      read (lun,*) iagro_thmb
      read (lun,*) agro_thmb_path
      read (lun,*) commu_MaxWait
      read (lun,*) commu_sleep
      read (lun,*) ffact
      read (lun,*) management_prescribed
      read (lun,*) planting_input
      read (lun,*) cultivar_input
      read (lun,*) overveg
      read (lun,*) isoilay
      read (lun,*) co2init
      read (lun,*) o2init
      read (lun,*) dtime
      read (lun,*) idiag
      read (lun,*,end=99) snorth
      read (lun,*,end=99) ssouth
      read (lun,*,end=99) swest
      read (lun,*,end=99) seast
      read (lun,*,end=99) flg_udg
      read (lun,*,end=99) udg_is,udg_ie
      read (lun,*,end=99) udg_js,udg_je
      read (lun,*,end=99) lon_llcorner,lat_llcorner
      read (lun,*,end=99) flg_debug, debug_ncounts
      read (lun,*,end=99) flg_detail
      read (lun,*,end=99) flg_dailyin_freq
      read (lun,*,end=99) flg_dailyout_freq
      read (lun,*,end=99) flg_daily4thmb
      read (lun,*,end=99) flg_sandlayers
      read (lun,*,end=99) flg_cmask
      read (lun,*,end=99) fileName_cmask
      read (lun,*,end=99) varName_cmask
      read (lun,*,end=99) flg_coverCrops
      read (lun,*,end=99) flg_mngt
      read (lun,*,end=99) myprocDir
c
 99   close (lun)
c
!-- pre-proc necessary information for the simulation

      irestart_glo = irestart
      irstyear_glo = irstyear
      iyear0_glo   = iyear0
      nrun_glo     = nrun
      iyearout_glo  = iyearout
      imonthout_glo = imonthout
      idailyout_glo = idailyout

      !make sure <agro_thmb_path> has '/' ended
      length = len(trim(agro_thmb_path))
      if( agro_thmb_path(length:length) .ne. '/') then
         agro_thmb_path(length+1:length+1) = '/'
      end if

      !make sure myprocDir only activate for mpi version
      if(nargc <=1 ) then  !not read-in user specified ibis.infile,this is a single processor run
        write(6,'(\,A)') '** This is a single processor run (because no specification of the optional commandline argument -f)'
        write(6,*) 'Ignoring the readin variable <myprocDir> from ibis.infile: ',trim(myprocDir)
        write(6,'(A,\)') 'Changing it back to ./'

        myprocDir = './'

      else

        flg_mpi = 1  !mpi run
        !get the complete simulation domain info from <proc.info>
        open(9999, status = 'old', file='proc.info')
          read(9999,*) str_buf,str_buf,istart_all,iend_all
          read(9999,*) str_buf,str_buf,jstart_all,jend_all
        close(9999)

      end if

      !ghost_io
      if(flg_ghost_io == 1) then !prefix "/tmp" to myproc for output directories (see command line processing).
        write(6,'(/,A)') 'ghost_io is requested from commandline'

        if(idailyout == 1 .or. ihourlyout == 1) then
          write(6,'(A)') 'Error! current version does not support ghost_io with <idailyout> or <ihourlyout> == 1'
          write(6,'(A)') 'make sure flg_dailyout_freq is enabled if daily data is needed'
          write(6,*) 'Stopping the program...'
          stop
        end if

        myprocDir_ori = myprocDir
        myprocDir = '/tmp/'//myprocDir  !the only place needed to specify /tmp or any desired temporary path

        call execute_command_line('mkdir -p '//trim(myprocDir)//'restart', 
     *                             EXITSTAT = exit_s,CMDSTAT = cmd_s, CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 'when creating ghost_io restart directory')

        call execute_command_line('mkdir -p '//trim(myprocDir)//'/output/yearly',
     *                             EXITSTAT = exit_s,CMDSTAT = cmd_s, CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 'when creating ghost_io yearly directory')

        call execute_command_line('mkdir -p '//trim(myprocDir)//'/output/monthly',
     *                             EXITSTAT = exit_s,CMDSTAT = cmd_s, CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 'when creating ghost_io monthly directory')

        call execute_command_line('mkdir -p '//trim(myprocDir)//'/output/daily',
     *                             EXITSTAT = exit_s,CMDSTAT = cmd_s, CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 'when creating ghost_io daily directory')

        call execute_command_line('mkdir -p '//trim(myprocDir)//'/output/hourly',
     *                             EXITSTAT = exit_s,CMDSTAT = cmd_s, CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 'when creating ghost_io hourly directory')

        write(6,'(A)') 'prefix /tmp to the output directories, create the corresponding directories complete...'
          
      end if

c tell user about the simulation
c
      write (*,*) ' '
      write (*,*) '*********************************************'
      write (*,*) '* Agro-IBIS: Integrated BIosphere Simulator *'
      write (*,*) '* IBIS U.S./Agroecosystem  Version          *'
      write (*,*) '*                                           *'
      write (*,*) '* V1.0 - Chris Kucharik                     *'
      write (*,*) '* kucharik@wisc.edu                         *'
      write (*,*) '* Center for Sustainability and             *'
      write (*,*) '* the Global Environment                    *'
      write (*,*) '*                                           *'
      write (*,*) '* University of Wisconsin-Madison           *'
      write (*,*) '*                                           *'
      write (*,*) '* Sept 1, 2005                              *'
      write (*,*) '*********************************************'
      write (*,*) ' '
c
      call date_and_time(p_date,p_time,p_zone,p_values) !Beginning of simulation time info
      write(6,'(/,8x,A,A,A,3(i2.2,A),/)') 'Program begins on ',p_date,' at local time ',p_values(5),':',p_values(6),':',p_values(7)

      !pre-proc management
      call preproc_management()

      if(flg_cmask == 1) then
         write(6,'(/,1x,A)') 'Using cmask file: '//trim(fileName_cmask)//'; variable name: '//trim(varName_cmask)
         if(flg_udg == 0) then
           write(6,*) 'Warning! not using udg. Simulation may fail (non UDG scenario not tested)'
         end if !on if flg_udg
      end if !on if flg_cmask

      if(flg_coverCrops == 1) then
         write(6,'(/,1x,A)') 'Enabling cover crops after corn/soy rotation harvest'
      end if !on if flg_coverCrops

      if(isoil_layer .eq. 1) then
         write(6,'(/,1x,A)') 'Enabling Read-in soils with layered data from a NetCDF file'
      end if

      if(iyr_pft .eq. 1) then
         write(6,'(/,1x,A)') 'Enabling Read-in land use on a year-by-year basis from NetCDF files'
      end if

      if(iyr_fert .eq. 1) then
         write(6,'(/,1x,A)') 'Enabling Read-in fertilizer application on a year-by-year basis from NetCDF files'
      end if

      if(iagro_thmb .eq. 1 ) then
         write(6,'(/,1x,A)') 'Enabling coupled Agro-IBIS THMB simulation'
         write(6,'(/,1x,A)') 'Mode: 1-way coupled (Agro-IBIS send data to THMB, but does not take THMB input)'
         write(6,'(1x,A,A)') 'Communication path at: ', trim(agro_thmb_path)
      end if

      if(iagro_thmb .eq. 2 ) then
         write(6,'(/,1x,A)') 'Enabling coupled Agro-IBIS THMB simulation'
         write(6,'(/,1x,A)') 'Mode: 2-way coupled (Agro-IBIS send data to THMB, and take THMB input)'
         write(6,'(1x,A,A)') 'Communication path at: ', trim(agro_thmb_path)
      end if

      if(flg_sandlayers == 1) then
         write(6,'(/,1x,A)') 'flg_sandlayer == 1: forcing last 2 layers as sand'
      end if 

      if (irestart .eq. 1) then
        write (*,*) 'running in restart mode'
      end if
c
      write (*,*) ' '
      write (*,9000) nrun
      if (.not. weather_generator) then
        write (*,*) 
     >  "ignoring iyranom, since we're not using the weather generator"
      else
        write (*,9005) iyranom
      end if
      write (*,*)
      write (*,*) 'Prescribed management? ', management_prescribed
      if (management_prescribed) then
        write (*,*) 'Planting_input: ', trim(planting_input)
        write (*,*) 'Cultivar_input: ', trim(cultivar_input)
      end if
      write (*,*) ' '

      write (*,9010) xres, yres
      write (*,9020) nlon, nlat

      if(flg_udg == 1 .or. flg_udg == 2) then
        write (*,*) ' '
        write(6,*) '**Running with user-defined grid**'
        write(6,*) 'flg_udg: ',flg_udg
        if(flg_udg == 1) then
          write(6,*) 'UDG: using lower left corner base'
        elseif(flg_udg == 2) then
          write(6,*) 'UDG: using upper left corner base'
        end if
        write(6,*) 'base corner (lon,lat):  ',lon_llcorner,lat_llcorner
        write(6,*) 'sub-grid longitudinal start and end index: ',udg_is,udg_ie
        write(6,*) 'sub-grid latitude     start and end index: ',udg_js,udg_je

        !simple check for grid number
        if( (udg_ie-udg_is+1)*(udg_je-udg_js+1) .ne. npoi) then
          write(6,*) 'Error! expecting number of grid points specified in ibis.infile equal <npoi> in <comgrid.f>'
          write(6,*) '(udg_ie-udg_is+1)*(udg_je-udg_js+1) = ', (udg_ie-udg_is+1)*(udg_je-udg_js+1)
          write(6,*) 'npoi = ', npoi
          write(6,*) 'Stopping the program...'
          stop
        end if

      end if

      write (*,*) ' '
c
 9000 format (1x,'length of this simulation (years)   : ',i8)
 9005 format (1x,'year to begin using anomalies       : ',i8)
 9006 format (1x,'number of iterations per day        : ',i8)
 9010 format (1x,'model lon, lat resolution (degrees) : ',2f8.5)
 9020 format (1x,'model domain (nlon x nlat)          : ',i5,
     >        ' by ',i5)  
 9030 format (1x,'last year run in this sequence      : ',i8)
c
c ---------------------------------------------------------------
c                     t i m e   c o n t r o l 
c ---------------------------------------------------------------
c
c determine the number of timesteps per day
c
      niter = int (86400.0 / dtime)
c
c test on length of the time step and total number of iteration per day
c
      test = amod (86400.0, dtime)
c
      if (test .gt. 1.e-20) then
        write (*,*) 'dtime ', dtime, ' should be divisible into 86400'  
        stop0
      else
        write (*,9006) niter
      end if
c
c check for restart conditions and set "iyrlast" accordingly, open
c ibis.annual in append mode
c
      if (irestart .eq. 1) then
        open (13, status='old', file='yearsrun.dat', err=9050)
        read (13,*) iyrlast
        close (13)
        goto 9059
 9050   write (*,*) 'warning: ibis.infile restart flag=1, indicating'
     >   //' restart, but yearsrun.dat not found'
        write (*,*) 'beginning from ', iyear0
        iyrlast = iyear0 - 1
 9059   open(20,file=trim(myprocDir)//'ibis.out.global',status='unknown',
     >       position='append')
        open(30,file=trim(myprocDir)//'ibis.out.vegtype',status='unknown',
     >       position='append')
      else
        iyrlast = iyear0 - 1
      end if
c
 
      iyrlast_glo = iyrlast  !for global communication, applied only once at start

      write (*,9030) nrun + iyrlast
      write (*,*) ' '
c
c ---------------------------------------------------------------
c                   i n i t i a l i z a t i o n
c ---------------------------------------------------------------
c
c calls to read in parameter files
c
      call rd_param 
      call read_input_description(datasource)

      !nc file info
      write(6,'(A,i8)') '(input_description <nc_file_bufrsize_KB>) Buffer size (KB) for nc file: ',nc_file_bufrsize_KB
c
c do some user input checks
c
      call check_inputs(nerr)
      if (nerr .gt. 0) then
        write(*,*)
        write(*,*) nerr, 'total errors found: HALTING'
        stop 1
      end if
      
c
c initial concentration of co2, o2
c
      co2conc = co2init
      o2conc  = o2init
c
c read global boundary condition datasets
c
      call readit (isimveg,snorth,ssouth,swest,seast,iwest,jnorth)
c
c allocate dynamically-allocated variables
c
      call allocate_vars(niter)
c
c check if diagnostic output is requested, if so read info from 'diag.infile'
c
      if (idiag .ne. 0) call inidiag(idiag)
c
c preliminary analysis of climate data
c
      if (irestart .eq. 0) then
        call climanl
      end if
c
c if we're using the weather generator, and if first year of this
c run/restart is an anomaly year, then read in anomalies for month 12 of
c previous year and month 1 of this year (note: in these calls, rdanom
c reads imonth+1, since next_month=TRUE) (note: we only do this if we're
c using the weather generator, because that is what needs the previous &
c next month of climate data)
c
      iy2 = iyrlast + nrun
      if (weather_generator .and. .not. read_daily_directly .and.
     >     (iyrlast+1 .ge. iyranom)) then
        call rdanom (11, iyrlast, iy2, iwest, jnorth, .TRUE.)
        call rdanom (12, iyrlast, iy2, iwest, jnorth, .TRUE.)
      end if
c
c 
      cropsums = 0
      if(iyr_pft == 1) then
        call initial_pft(irestart,iyear0,irstyear)  ! initialize pft including imaize, isoybean, iwheat
      end if

      cropsums = imaize + isoybean + iwheat + irotation
      if (iwheat .gt. 0) then
        iwheat_index = iwheat
      else  ! iwheat .le. 0, in which case iwheat_index doesn't matter
        iwheat_index = 1
      end if
c
c initialize the model
c
      call initial (isimveg, irestart, iyrlast)
c
c initialize crop variables
c have to be initialized when crops are replacing natural vegetation in
c a series of simulations
c
      if (cropsums .gt. 0) then
        call initialcrop
        call read_cfert
      end if
  
c
c initialize random number generator seed
c
      if (seed .eq. 0) then
        seed = random_seed()
      end if
      write (*,*) ' '
      write (*,*) 'Seeding the random number generator with: ', seed
      write (*,*) ' '
      dran = ran2 (seed)

c initialize io-dyn variables

      if(flg_dailyin_freq == 1 .or. flg_dailyin_freq == 2) then
        call alloc_dailyin_variables
      end if

      if(flg_dailyout_freq == 1 .or. flg_dailyout_freq == 2) then
        call alloc_dailyout_variables
        if(flg_daily4thmb == 1 .or.flg_daily4thmb == 2) then
          call alloc_daily4thmb_variables
        end if !on if flg_daily4thmb
      end if !on if flg_dailyout_freq

c
c ---------------------------------------------------------------
c              s t a r t   o f   t i m e   l o o p s
c ---------------------------------------------------------------
c
      write (*,*) ' '
      write (*,*) '********************************'
      write (*,*) '* start of the IBIS simulation *'
      write (*,*) '********************************'
      write (*,*) ' '
c
c reset elapsed number of days, accumulating past days if restart mode
c
      nday = 0
c
      if (irestart .eq. 1) then
        do 150 iyear = iyear0, iyrlast
          nday = nday + 365
          if (is_leap(iyear)) then
            nday = nday + 1
          end if
 150    continue
      end if
      ndayr=nday
c
c start of yearly loop
c
      iy1 = iyrlast + 1
      iy2 = iyrlast + nrun

c
      do 200 iyear = iy1, iy2

        call system_clock(yearly_ts,clock_rate,clock_max)

        year_idx = iyear - iyrlast
c
c reset julian date
c
        jday = 0

        if(flg_dailyout_freq == 2) dyn_idx = 0 !reset index

        if((iyear.ge.iyrdaily .or. .not. weather_generator) .and. 
     *       (iyear .lt. imetyear .or. iyear .gt. imetend)) then
          flg_read_daily = .true.
        end if

c
c determine the calendar for this year: leap vs. non-leap year
c
        if (is_leap(iyear)) then
          ndaypm(2) = 29
          ndaypy = 366
        else
          ndaypm(2) = 28
          ndaypy = 365
        end if          
c
c if we're not using the weather generator, and this year is prior to
c iyrdaily, then choose a random year of climate data to use this year
c
        if (iyear .lt. iyrdaily .and. .not. weather_generator) then
          if (is_leap(iyear)) then
            call random_choice(iyear_climate, avail_leaps, size(avail_leaps), seed)
          else
            call random_choice(iyear_climate, avail_nonleaps, size(avail_nonleaps), seed)
          end if
          write(*,*)
          write(*,'(i4, ": Using ", i4, " from climate data")') iyear, iyear_climate

        else  ! either iyear >= iyrdaily, or we're using the weather generator
          iyear_climate = iyear
        end if

        if(flg_mngt == 1 .and. (.not. climate_mngt_mode=='none')) then !climate management
          iyear_climate = climateMngt_yearsList(year_idx)
          write(*,'(i4, ": climate management using ", i4, " from climate data")') iyear, iyear_climate
        end if !on if flg_mngt & climate_mngt_mode

c
c JCP (8-17-10): Open hourly input/output files for this year if using hourly met data
c
        if ((iyear .ge. imetyear) .and. (iyear .le. imetend)) then
           call openHourlyMet(iyear, seed, soilcspin)
           call openHourlyOut
        endif
c
c Y.Li (2014-06-23): added hourly output when ihourlyout = 1
c
        if ( ihourlyout .eq.1 ) then

          if ((iyear .ge. imetyear) .and. (iyear .le. imetend)) then
            !in case both hourly met data and ihourlyout = 1 are enabled
            call closeHourlyOut
          else
            call openHourlyOut
          end if

        end if !on if ( ihourlyout .eq.1 )
c
c read this year's management data
c
        if (management_prescribed) then
          call read_management_data('input', iyear, iwest, jnorth, planting_input, cultivar_input)
        end if

c read nc when dailyout_freq enabled & read daily is needed
        if(flg_read_daily == 1 .and. flg_dailyout_freq == 2) then !chunk data for the year
          call readnc_dailyout_freq(imonth,iyear,iwest,jnorth)
        end if
c
c start of monthly loop
c 
         do 210 imonth = 1, 12
c
          if (.not. read_daily_directly .and.
     >          (iyear .lt. imetyear .or. iyear .gt. imetend)) then
            if (.not. weather_generator) then
c When weather_generator=FALSE, we only need the current month's data
c (and it would be messy to try to read the previous & next month's data
c when we're using a random year of climate data); so here we simply
c read the current month's data. Also, with weather_generator=FALSE, we
c ignore the value of iyranom: we will ALWAYS read anomalies for the
c year given by iyear_climate.
              call rdanom (imonth, iyear_climate, iy2, iwest, jnorth, .FALSE.)
            else if (iyear_climate .ge. iyranom .or.
     >        (iyear_climate .eq. iyranom-1 .and. imonth .ge. 11)) then
c When weather_generator=TRUE, we need the previous and next month's
c data as well as the current month's data; so here we read the next
c month's data, and start reading anomalies with the December before
c iyranom. For simplicity, we always read next month's data when
c weather_generator=TRUE, even if iyear >= iyrdaily so that we're no
c longer using the weather generator. Note also that in this case
c iyear_climate = iyear.
              call rdanom (imonth, iyear_climate, iy2, iwest, jnorth, .TRUE.)
            end if
          endif

c read nc when dailyout_freq enabled
         if(flg_read_daily == 1 .and. flg_dailyout_freq == 1) then !chunk data for the month
           call readnc_dailyout_freq(imonth,iyear,iwest,jnorth)
         end if

         if(flg_dailyout_freq == 1) dyn_idx = 0  !reset index

c
c start of daily loop
c
         do 220 iday = 1, ndaypm(imonth)


           call system_clock(clock_s,clock_rate,clock_max)
c
c update user on model progress once a day
c
            if(iday == 1) then
              write (*,9100) iyear, imonth, iday
            else
              write (*,9100,advance='no') iyear, imonth, iday
            end if
 9100       format (1x,'starting - year/month/day: ',i4,'/',i2,'/',i2)
c
c update elapsed number of days and julian date
c
            nday = nday + 1
            jday = jday + 1

            dyn_idx = dyn_idx + 1
c
c pre-calculate coszen for each latitude and time, 
c and the daily-integrated coszen for each latitude
c
            call calc_coszen(jday)
c
c get daily means (if requested)
c and determine today's climatic conditions
c
            if ((iyear.ge.iyrdaily .or. .not. weather_generator)
     >        .and. (iyear .lt. imetyear .or. iyear .gt. imetend)) then
c don't use weather generator: read today's climate data
              if(flg_read_daily == 1 .and. (flg_dailyout_freq == 1 .or. flg_dailyout_freq == 2) ) then  
                call assign_dailyin_freq  !assign read-in chunk data to daily variables
              else
                call rdday(jday,imonth,iyear_climate,iwest,jnorth)
              end if !on if flg_dailyout_freq

              call daily (iyear, imonth, iday, jday, seed, 1, iyrlast, nrun)
            else if (iyear   .ge. imetyear .and.
     >               iyear   .le. imetend  .and.
     >               iday    .eq. 1        .and.
     >               imonth  .eq. 1        ) then
              call dailymet (imonth, iday, seed, 1)
            else if (iyear .lt. imetyear .or. iyear .gt. imetend) then  
c use weather generator
              call daily (iyear, imonth, iday, jday, seed, 0, iyrlast, nrun)
            end if
c
c determine the daily vegetation cover characteristics
c
            call pheno(jday)

c modified for croplands
c 1-17-02 CJK
c if rotation of crops is occuring, call on first day of year
c
c added read-in pft data in a yearly basis, Y.Li 2014-06-18
            if ( (irotation .gt. 0 .and. jday .eq. 1).or.
     >           (iyr_pft .eq. 1   .and. jday .eq.1) ) 
     >          call rotation(irestart, irstyear, iyear)
c
c determine if planting is able to take place for crop and location
c
      cropsums = 0
      cropsums = imaize + isoybean + iwheat + irotation 
c
      if ( (cropsums .gt. 0).or.
     >     (iyr_fert.eq.1  )   )

     >    call planting(irestart, irstyear, iyear0,iyear,imonth,iday, jday,ffact)
c
c call daily crop phenology subroutine
c
      if (cropsums .gt. 0)
     >    call phenocrop(iyear, iyear0, imonth, iday, jday)
c
c call soil biogeochemistry model
c
c soil spin up procedure calls soil bgc model spincons times
c for each time step to spin up the process (accelerate C & N pools)
c
            spinfrac  = 0.75          ! fraction of nspinsoil time used to
c                                     ! spin up soil at maximum spin up rate
c
            spincons  = 40.0          ! number of times soilbgc subroutine is
c                                     ! called for each day in simulation during
c                                     ! acceleration procedure
c
            eqyears   = 50            ! years used past spin up to slowly bring
c                                     ! entire soil to a steady equilibrium
c
            nspinsoil = iyear0 + 150  ! year of simulation that soil c reaches equilibrium
c
c
            if (soilcspin .eq. 1) then

              if ((iyear - iyear0) .le.
     >            (spinfrac * (nspinsoil - iyear0 - eqyears))) then
                 spinmax = int(spincons)
c
              else if ((iyear - iyear0) .lt.
     >               (nspinsoil - iyear0 -  eqyears)) then
c
                slope   = spincons / ((nspinsoil - iyear0 - eqyears) -
     >                    (spinfrac * (nspinsoil - iyear0 - eqyears)))
c
                spinmax = int (spincons - (slope * ((iyear - iyear0) -
     >            (spinfrac * (nspinsoil - iyear0 - eqyears)))))
c
                spinmax = max(spinmax,1)
c
              else
c
                spinmax = 1
c
              endif
c
            else 
c
                spinmax = 1
c
            endif
c
            do 230 spin = 1, spinmax
              call soilbgc (iyear, iyear0, imonth, iday, jday, nspinsoil,
     >                      spin, spinmax)
 230        continue
c
c determine the length of a precipitation event (between 4 and 24 hours),
c and time to start and end precipitation event. plen is in timesteps, while
c plens, startp, and endp are in seconds
c
            plenmin = 1 +  int ((4.0 * 3600. - 1.) / dtime)
            plenmax = max (int (24.0 * 3600. / dtime), plenmin)
c
            if (iyear .lt. imetyear .or. iyear .gt. imetend) then
              plen    = min (plenmax, int (
     >                     plenmin + ran2(seed) * (plenmax-plenmin+1) ))
cjk
            else
              plen = 1
            endif
c
c plen = 1 for hourly precipitation events
c 
            plens   = dtime * plen
            startp  = dtime * min (niter-plen, 
     >                             int(ran2(seed)*(niter-plen+1)))
            endp    = startp + plens
c
c calculate irrigation amount, timing, and duration
c if applicable to model run
c
c assume irrigation duration for a day is 12 hours long between 
c 6 am and 6 pm local time - this might be changed if managed irrigation can
c take place at night (e.g. crops)  
c
             ilens  = dtime * (12.0 * 3600. / dtime)   
             starti = dtime * (6.0  * 3600. / dtime)
             endi   = starti + ilens 

             call irrigation(iday,imonth) 

             if(iagro_thmb == 1 .or. iagro_thmb == 2) then ! Y.Li (enable daily coupled simulation for agro-ibis and thmb)

               call agro_thmb_read_input(ndaypm,iyear,iyear-iy1,imonth,iday,thmb_lakem,'lakem','lakem.wat.nc')

             end if

c
c start of hourly loop
c
            do 240 istep = 1, niter
c
c calculate the time of day since midnight (in seconds)
c
              time = get_time(istep)

              it_hour = istep
c
c determine climatic conditions for the given timestep
c also apply irrigated water
c
c check to see if using hourly met input data, otherwise run diurnal sub
c to determine conditions for timestep
c
              if (iyear .ge. imetyear .and. iyear .le. imetend) then

c JCP (8-17-10): Commenting out old hourly met file stuff
c                call inimet (iyear,imonth,iday,nday,time,seed)
c                call diurnalmet(istep, jday, plens, startp, endp, seed,
c     >                      ilens, starti, endi)
c
c JCP (8-17-10): Adding in my own hourly met file stuff
c
                 call readHourlyMet(iyear, imonth, iday, istep,
     >                jday, irrigate_flag, ilens, starti, endi)

              else
                call diurnal (istep, jday, plens, startp, endp, seed,
     >                      ilens, starti, endi)
              endif
c
c -------------------------------------------------------------------
c added for Green-Ampt infiltration model
c              t_sec = time
c              t_startp = startp
c -------------------------------------------------------------------
c call the land surface model
c
              call lsxmain
c
c accumulate some variables every timestep
c
              call sumnow
c
              call sumday   (istep,plens,iyear,imetyear,imetend) 
              call summonth (istep, iday, imonth)
              call sumyear  (istep, iday, imonth)
c
c reset length of fire season at beginning of every year, 
c otherwise keep running total- MMM 4/2010
c
              do i = 1, npoi  
                if(cmask(i) == 0 ) cycle

                 if (jday.eq.1) firelength(i) = 0
                 if (istep.eq.niter) then
                    firelength(i) = firelength(i) + adprobfire(i)
                 endif
              enddo
              
c
c call to nitrogen stress routine for crops
c 
         call nitrostress(istep, iday,imonth)
c       
         call leaching(irestart, irstyear, istep,iday,imonth,iyear,
     >                 iholdsoiln, iyear0)

c
c JCP (8-17-10): Write hourly output if using hourly input data
c
c Y.Li (2014-06-23) or when ihourlyout = 1
c
         if ( (iyear .ge. imetyear .and. iyear .le. imetend).or.
     >        (ihourlyout.eq.1) ) then
            call whourly(iyear, imonth, iday, istep)
         endif

        if(flg_dailyout_freq == 1 .or. flg_dailyout_freq == 2) then
          call get_dailyout_freq(istep,niter)
        end if

c
c write out diagnostics
c
              if (idiag .ne. 0) then
c
                do 250 i = 1, idiag
                  if (iyear.ge.diagstart(i) .and.
     >                iyear.le.diagend(i)   .and. 
     >                ndiagpt(i).ge.1)  then  
                    if (mod(istep,nfreq(i)).eq.0) then
                      call wdiag (i, iyear, imonth, iday, istep)
                    end if
                  end if
 250            continue
c  
              end if
c  
c write out hourly canopy photosynthesis rates
c 
c           open(33,file='canopy.ps.dat',status='unknown')
c           if (iyear .eq. 1998 .and. (imonth .gt. 3 .and. imonth .lt. 10)) then
c              write(33,270) iyear, imonth, iday, istep, ancub(4)*1.e+06,
c     >        ancuc(4)*1.e+06, totlaiu(4), plai(4,5), lai(4,2)
c              write(33,270) iyear, imonth, iday, istep, agcc4(4)*1.e+06,
c     >        ancc4(4)*1.e+06, plai(4,4), plai(4,5), plai(4,11), plai(4,12)
c               write(33,270) iyear, imonth, iday, istep, ancc4(4)*1.e+06*totlail(4),
c     >                       lai(4,1), frac(4,14), fl(4) 
c            endif
c 270        format(i5,i3,i3,i4,4f7.2)
              
c
c ---------------------------------------------------------------
c               e n d   o f   t i m e   l o o p s
c ---------------------------------------------------------------
c
c end of the hourly loop
c
 240        continue
c
c write out daily output
c
            if (idailyout.eq.1) then
c             write (*,*) ' '
c             write (*,*) 'writing daily output'
c             write (*,*) ' '
              call wdaily (nday, iyear, iyear0, ndayr, jday, irstyear, irestart) 
            endif

            if(flg_daily4thmb == 1 .and. flg_dailyout_freq == 0) then

              call  write_dailyOutput4THMB(jday,iyear,imonth,iday)

            end if !on if(flg_daily4thmb)

            if(iagro_thmb == 1 .or. iagro_thmb == 2) then

              call agro_thmb_write_output(iyear,imonth,iday,adsrunoff,'adsrunoff','adsrunoff.nc',0)    !daily average surface runoff
              call agro_thmb_write_output(iyear,imonth,iday,addrainage,'addrainage','addrainage.nc',0) !daily average drainage
              call agro_thmb_write_output(iyear,imonth,iday,adrain,'adrain','adrain.nc',0)             !daily average rain
              call agro_thmb_write_output(iyear,imonth,iday,adevap,'adevap','adevap.nc',1)             !daily average evaporation

            end if


          !report performance of a day
          call system_clock(clock_e,clock_rate,clock_max)
          write(6,'(A,g16.8,A)') ' **run time of the day: ',real(clock_e-clock_s)/(real(clock_rate)),' seconds'

          !end the simulation for debug mode irestart = -1 (daily)
          if(flg_debug == 1 .and. iday <= debug_ncounts) then
            write(6,*)
            write(6,*) 'Warning! This is a debug mode flg_debug = 1 (daily) at day#: ',iday
            write(6,*) 'Put the simulation to end after n-day simulation after reporting results: ',debug_ncounts

            call gdiag (iyear, iyear0)
            call vdiag (iyear, iyear0)

            if(iday == debug_ncounts) then
              goto 200   !end of the yearly loop
            end if
          end if  !on if(flg_debug)
          !end debug mode

c
c end of the daily loop
c
 220      continue
c
c write out monthly output (use 1st day of this month because of GrADS bug)
c
          if (imonthout.eq.1 .or. imonthout.eq.2) then
            write (*,*) ' '
            write (*,*) 'writing monthly output'
            write (*,*) ' '
            call wmonthly (nday-ndaypm(imonth)+1, imonth, iyear, 
     >                     iyear0, irstyear, jday, irestart)
          endif

          !end the simulation for debug mode irestart = -2 (month)
          if(flg_debug == 2 .and. imonth <= debug_ncounts) then
            write(6,*)
            write(6,*) 'Warning! This is a debug mode flg_debug = 2 (month) at month#: ',imonth
            write(6,*) 'Put the simulation to end after n-month simulation after reporting results: ',debug_ncounts

            call gdiag (iyear, iyear0)
            call vdiag (iyear, iyear0)
            if(imonth == debug_ncounts) then
              goto 200   !end of the yearly loop
            end if !on if(imonth)
          end if !on if (flg_debug)
          !end debug mode

          !write daily output in a monthly frequency if requested
          if(flg_dailyout_freq == 1) then
            call write_dailyout_dyn(imonth,iyear)
          end if

c
c end of the monthly loop
c
 210    continue
c
c recalculate bioclimatic parameters if using anomalies
c Note: if read_daily_directly, then we will call this unnecessarily for iyear = iyranom-1,
c but that's not a problem.
c
        if (.not. weather_generator .or.
     >    (iyear .ge. iyranom-1 .or. iyear.ge.iyrdaily)) then
          call climanl2
        end if
c
c get new o2 and co2 concentrations for this year
c
        if (isimco2.eq.1) call co2 (co2init, co2conc, iyear)
c
c perform vegetation dynamics
c
        if (isimveg.ne.0) call dynaveg (iyear, isimfire)
c
c calculate simple annual diagnostic variables for the globe
c
        call gdiag (iyear, iyear0)
c
c calculate simple annual diagnostic variables for each vegetation type
c
        call vdiag (iyear, iyear0)
c
c write out annual output
c
        if (iyearout.eq.1 .or. iyearout.eq.2) then
c         write (*,*) ' '
c         write (*,*) 'writing annual output'
c         write (*,*) ' '
          call wyearly (nday, iyear, iyear0, irstyear, irestart)
        endif

        !write daily output in a yearly frequency if requested
        if(flg_dailyout_freq == 2) then
          call write_dailyout_dyn(imonth,iyear)
        end if
c
c write restart files
c
! (Y.Li)
! If you don't want to write restart files, e.g., debuging using the
! same case, or trying different methodologies for the same period of
! time, and you want to avoid copying the restart files and changing
! yearsrun.dat --> set flg_wrestart = 0 in ibis.infile

        if (flg_wrestart == 1) then

          call wrestart (nday, iyear, iyear0)
c
c update iyrlast value and file yearsrun.dat
c
          iyrlast = iyrlast + 1
          open(12,file=trim(myprocDir)//'yearsrun.dat',status='unknown')
          write(12,*) iyrlast
          close(12)

        end if

c JCP (8-17-10): Close hourly data input/output files
        if ((iyear .ge. imetyear) .and. (iyear .le. imetend)) then
           call closeHourlyMet()
           call closeHourlyOut
        endif
c
c Y.Li (2014-06-23): close hourly output when ihourlyout = 1
c
        if ( (ihourlyout .eq.1) .and.
     >       (.not.((iyear .ge. imetyear) .and. (iyear .le. imetend))) )then
          call closeHourlyOut
        end if

        if( flg_ghost_io == 1) then  !tar and move daily and hourly data back to normal output directories
          write(num,'(i4.4)')  iyear

          ! Start of changes by Cicada Dennis April-June 2016.
          ! Changed the moving of data so it uses a tar pipe, rather than writing tar files to disk.
          ! Also added a flag to add the writing of diagnostic information on the disk usage of the
          ! node's ramdisk - /tmp (ghost_io) file system.
          if ( flg_sys_stats == 1) then
            !print out disk usage now (added to diagnose full file system)
            write(0,*) 'About to list dirs before transfer, year: ', iyear
c           The commented out code is implemented in fs_report.sh
c           This was done in order to stop errors that occasionally occur
c           when getting the disk usage, from terminating the ibis process.
c            call execute_command_line('ls -l '//trim(myprocDir)//
c     *                                'output/* >&2 && df >&2 && '//
c     *                                'du -d 0 /tmp/mpi >&2 '//
c     *                                ' && du '//trim(myprocDir)//' >&2',
c     *                                 EXITSTAT = exit_s,
c     *                                 CMDSTAT = cmd_s, 
c     *                                 CMDMSG = message)
            call execute_command_line('echo "Working directory:" >&2'//
     *                                ' && pwd >&2 && ./fs_report.sh '//
     *                                trim(myprocDir)//'output/\* '//
     *                                '/tmp/mpi '//
     *                                trim(myprocDir)//' >&2',
     *                                 EXITSTAT = exit_s,
     *                                 CMDSTAT = cmd_s, 
     *                                 CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *                 'when trying to print directory listings before'//
     *                 ' deleting files in: '//trim(myprocDir)//
     *                 ' at year:'//num)
          end if !on if flg_sys_stats

          !tar the daily files
          write(0,*) 'About to tar daily files.'
          call execute_command_line(' (cd '//trim(myprocDir)//
     *                              'output/daily/ && '//
     *                              'tar -c --ignore-failed-read *) |'//
     *                              ' (cd '//trim(myprocDir_ori)//
     *                              'output/daily/ && tar -xf -)', 
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
          call error_handler_command(exit_s,cmd_s,message, 
     *               'when moving ghost_io daily at year '//num)

          !remove the daily files from /tmp
          write(0,*) 'About to remove daily files.'
          call execute_command_line('rm -f '//trim(myprocDir)//
     *                              'output/daily/*',
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
          call error_handler_command(exit_s,cmd_s,message, 
     *               'when deleting ghost_io daily at year '//num)

          if(flg_daily4thmb == 1 .or. flg_daily4thmb == 2) then
            !tar the hourly files
            write(0,*) 'About to tar hourly files.'
            call execute_command_line(' (cd '//trim(myprocDir)//
     *                              'output/hourly/ && '//
     *                              'tar -c --ignore-failed-read *) |'//
     *                              ' (cd '//
     *                              trim(myprocDir_ori)//
     *                              'output/hourly/ && tar -xf -)', 
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *                 'when moving ghost_io hourly at year '//num)

            !remove the hourly files from /tmp
            write(0,*) 'About to remove hourly files.'
            call execute_command_line('rm -f '//trim(myprocDir)//
     *                                'output/hourly/*',
     *                                 EXITSTAT = exit_s,
     *                                 CMDSTAT = cmd_s, 
     *                                 CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *                 'when deleting ghost_io hourly at year '//num)
          end if

          ! The transfer time of yearly and monthly files can be
          ! controlled through parameters within the ibis.infile.
          ! The transfer times can be controlled through iyearout and imonthout.
          ! When those parameters are set to 2, then the files are transfered 
          ! at the end of year instead of the end of the simulation.
          ! Restart files are not transfered until the end of the simulation.
          if (iyearout.eq.2) then
            !tar the yearly files on /tmp to dc2
            write(0,*) 'About to tar yearly files.'
            call execute_command_line(' (cd '//trim(myprocDir)//
     *                              'output/yearly/ && '//
     *                              'tar -c --ignore-failed-read *) |'//
     *                              ' (cd '//trim(myprocDir_ori)//
     *                              'output/yearly/ && tar -xf -)', 
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *               'when moving ghost_io yearly')

            !remove the yearly files from /tmp
            write(0,*) 'About to remove yearly files.'
            call execute_command_line('rm -f '//trim(myprocDir)//
     *                                'output/yearly/*',
     *                                 EXITSTAT = exit_s,
     *                                 CMDSTAT = cmd_s, 
     *                                 CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *               'when deleting ghost_io yearly at year '//num)
          end if ! (iyearout.eq.2)

          if (imonthout.eq.2) then
            !tar the monthly files on /tmp to dc2
            write(0,*) 'About to tar monthly files.'
            call execute_command_line(' (cd '//trim(myprocDir)//
     *                              'output/monthly/ && '//
     *                              'tar -c --ignore-failed-read *) |'//
     *                              ' (cd '//trim(myprocDir_ori)//
     *                              'output/monthly/ && tar -xf -)', 
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *               'when moving ghost_io monthly')

            !remove the monthly files from /tmp
            write(0,*) 'About to remove monthly files.'
            call execute_command_line('rm -f '//trim(myprocDir)//
     *                                'output/monthly/*',
     *                                 EXITSTAT = exit_s,
     *                                 CMDSTAT = cmd_s, 
     *                                 CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *               'when deleting ghost_io monthly at year '//num)
          end if ! (imonthout.eq.2)

          !tar the restart files on /tmp to dc2
          !This copied over them each time, but each time will have
          !files more recent than the last time copied.
          write(0,*) 'About to tar restart files '//
     *        '(they do not get removed, just copied over each time).'
          call execute_command_line(' (cd '//trim(myprocDir)//
     *                              'restart/ && '//
     *                              'tar -c --ignore-failed-read *) |'//
     *                              ' (cd '//trim(myprocDir_ori)//
     *                              'restart/ && tar -xf -)', 
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
          call error_handler_command(exit_s,cmd_s,message, 
     *             'when moving ghost_io restart')

          if ( flg_sys_stats == 1) then
            !print disk usage again to see what is left (added to diagnose full file system)
            write(0,*) 'About to list directories after data transfer.'
            call execute_command_line('./fs_report.sh '//
     *                                trim(myprocDir)//'output/\* '//
     *                                '/tmp/mpi '//
     *                                trim(myprocDir)//' >&2',
     *                                 EXITSTAT = exit_s,
     *                                 CMDSTAT = cmd_s, 
     *                                 CMDMSG = message)
            call error_handler_command(exit_s,cmd_s,message, 
     *                 'when trying to print directory listings after'//
     *                 ' deleting files in: '//trim(myprocDir)//
     *                 ' at year:'//num)
          ! End of changes by Cicada Dennis

          end if !on if flg_sys_stats

        end if !on if flg_ghost_io

c
c end of the yearly loop
c
        call system_clock(yearly_te,clock_rate,clock_max)
        write(6,'(A,g16.8,A)') ' **run time of the year: ',real(yearly_te-yearly_ts)/(real(clock_rate)),' seconds'

! Slavin - Close the input data files that we've kept open for this year.
        call close_logged_ncids
 
 200  continue
c
c end of the simulation
c
      if( flg_ghost_io == 1) then !deal with restart files
        !tar any remaining yearly files on /tmp to dc2
        write(0,*) 'About to tar yearly files at end of the simulation.'

        ! Start of changes by Cicada Dennis, April/May 2016
        call execute_command_line('cp '//trim(myprocDir)//
     *                            'yearsrun.dat '//
     *                            trim(myprocDir_ori),
     *                             EXITSTAT = exit_s,
     *                             CMDSTAT = cmd_s, 
     *                             CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 
     *               'when moving yearsrun.dat')

        call execute_command_line(' (cd '//trim(myprocDir)//
     *                            'output/yearly/ && '//
     *                            'tar -c --ignore-failed-read *) |'//
     *                            ' (cd '//trim(myprocDir_ori)//
     *                            'output/yearly/ && tar -xf -)', 
     *                             EXITSTAT = exit_s,
     *                             CMDSTAT = cmd_s, 
     *                             CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 
     *               'when moving ghost_io yearly')

        !tar any remaining monthly files on /tmp to dc2
        write(0,*) 'About to tar monthly files at end of simulation.'
        call execute_command_line(' (cd '//trim(myprocDir)//
     *                            'output/monthly/ && '//
     *                            'tar -c --ignore-failed-read *) |'//
     *                            ' (cd '//trim(myprocDir_ori)//
     *                            'output/monthly/ && tar -xf -)', 
     *                             EXITSTAT = exit_s,
     *                             CMDSTAT = cmd_s, 
     *                             CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 
     *               'when moving ghost_io monthly')

        !tar the files
        write(0,*) 'About to tar restart files at end of simulation.'
        call execute_command_line(' (cd '//trim(myprocDir)//
     *                            'restart/ && '//
     *                            'tar -c --ignore-failed-read *) |'//
     *                            ' (cd '//trim(myprocDir_ori)//
     *                            'restart/ && tar -xf -)', 
     *                             EXITSTAT = exit_s,
     *                             CMDSTAT = cmd_s, 
     *                             CMDMSG = message)
        call error_handler_command(exit_s,cmd_s,message, 
     *             'when moving ghost_io restart')

        if ( flg_sys_stats == 1) then
          !print disk usage again to see what is left (added to diagnose full file system)
          write(0,*) 'About to list dirs after simulation is over.'
          call execute_command_line('./fs_report.sh '//
     *                              trim(myprocDir)//'output/\* '//
     *                              '/tmp/mpi '//
     *                              trim(myprocDir)//' >&2',
     *                               EXITSTAT = exit_s,
     *                               CMDSTAT = cmd_s, 
     *                               CMDMSG = message)
          call error_handler_command(exit_s,cmd_s,message, 
     *               'when trying to print directory listings after'//
     *               ' deleting files in: '//trim(myprocDir)//
     *               ' at year:'//num)
        end if !on if flg_sys_stats
        ! End of changes by Cicada Dennis

      end if !on if flg_ghost_io

      !clean up 
      if(flg_dailyin_freq == 1 .or. flg_dailyin_freq == 2) then
        call dealloc_dailyin_variables
      end if

      if(flg_dailyout_freq == 1 .or. flg_dailyout_freq == 2) then
        call dealloc_dailyout_variables
        if(flg_daily4thmb == 1 .or. flg_daily4thmb == 2) then
          call dealloc_daily4thmb_variables
        end if !on if flg_daily4thmb
      end if !on if flg_dailyout_freq

      if(flg_mngt == 1) then
        call dealloc_management_variables
      end if !on if flg_mngt

      call system_clock(clock_te,clock_rate,clock_max)

      write(6,*) '*** total run time ***:',real(clock_te-clock_ts)/(real(clock_rate)),' seconds'

      write (*,*) ' '
      write (*,*) '*** end of run ***'
      write (*,*) ' '

      call date_and_time(p_date,p_time,p_zone,p_values) !End of simulation time info
      write(6,'(/,8x,A,A,A,3(i2.2,A),/)') 'Program ends on ',p_date,' at local time ',p_values(5),':',p_values(6),':',p_values(7)
c
      stop0

c ---------------------------------------------------------------
c Internal subroutines
c ---------------------------------------------------------------

      contains

c ---------------------------------------------------------------
      subroutine check_inputs(nerr)
c ---------------------------------------------------------------
c
c check user inputs, return # of errors found in nerr
c
c Subroutine arguments
c
      integer nerr

      nerr = 0

      if (seed .gt. 0) then
        write(*,*)
        write(*,*) 'ERROR: seed must be <= 0'
        write(*,*) '(where 0 signifies initialization using CPU time)'
        write(*,*) 'seed = ', seed
        nerr = nerr + 1
      end if

      if (iyrdaily .lt. istyrd) then
        write(*,*)
        write(*,*) 'ERROR: iyrdaily < istyrd'
        write(*,*) 'iyrdaily = ', iyrdaily
        write(*,*) 'istyrd = ', istyrd
        nerr = nerr + 1
      end if

      if (read_daily_directly) then
        if (iyranom .ne. iyrdaily) then
c     for read_daily_directly, we don't use iyranom; we make sure the
c     user knows this by ensuring that they set iyranom = iyrdaily
          write(*,*)
          write(*,*) 'ERROR: for read_daily_directly, you must set iyranom = iyrdaily'
          write(*,*) '(to signify that you know that iyranom is ignored)'
          nerr = nerr + 1
        end if 
      else  ! .not. read_daily_directly
        if (weather_generator) then
          if (iyranom .le. istyrm) then
            write(*,*)
            write(*,*) 'ERROR: for read_daily_directly=FALSE and weather_generator=TRUE,'
            write(*,*) 'iyranom must be > istyrm'
            write(*,*) '(strictly greater to allow reading of Dec. of previous year)'
            write(*,*) 'iyranom = ', iyranom
            write(*,*) 'istyrm = ', istyrm
            nerr = nerr + 1
          end if
        else  ! .not. weather_generator
          if (iyrdaily .lt. istyrm) then
            write(*,*)
            write(*,*) 'ERROR: for read_daily_directly=FALSE and weather_generator=FALSE,'
            write(*,*) 'iyrdaily must be >= istyrm (since we ignore iyranom in this case)'
            write(*,*) 'iyrdaily = ', iyrdaily
            write(*,*) 'istyrm = ', istyrm
            nerr = nerr + 1
          end if
        end if
      end if

      if (.not. weather_generator) then
        if (size(avail_leaps) .eq. 0 .or. size(avail_nonleaps) .eq. 0) then
          write(*,*)
          write(*,*) 'ERROR: for weather_generator=FALSE'
          write(*,*) '(i.e., using random yearsof climate data for spin-up),'
          write(*,*) 'there must be at least one available leap year'
          write(*,*) 'and one available non-leap year in the climate data set'
          write(*,*)
          write(*,*) '# avail_leaps = ', size(avail_leaps)
          write(*,*) '# avail_nonleaps = ', size(avail_nonleaps)
          nerr = nerr + 1
        end if
      end if

      return
      end subroutine check_inputs


      end
c
c
c ---------------------------------------------------------------
      subroutine lsxmain
c ---------------------------------------------------------------
c
c uses:
c 
      use comgrid
      use compar
      use com1d
      use comatm
      use comsoi
      use combcs,only: cmask
c
      implicit none
c
c Local variables
c
      integer ib,   ! waveband number (1= visible, 2= near-IR)
     >        i,    ! loop indice
     >        iday,
     >        imonth,
     >        iyear,
     >        iyear0
c
c ---------------------------------------------------------------
c added for Green-Ampt infiltration model
c      call initsw
c ---------------------------------------------------------------
c set physical soil quantities
c
      call setsoi
c
c calculate areal fractions wetted by intercepted h2o
c
      call fwetcal
c
c set up for solar calculations
c
      call solset
c
c solar calculations for each waveband
c
      do 100 ib = 1, nband
c
c solsur sets surface albedos for soil and snow
c solalb performs the albedo calculations
c solarf uses the unit-incident-flux results from solalb
c to obtain absorbed fluxes sol[u,s,l,g,i] and 
c incident pars sunp[u,l]
c
        call solsur (ib)
        call solalb (ib)
        call solarf (ib)
c
 100  continue
c
c calculate ir fluxes
c
      call irrad
c
c step intercepted h2o
c
      call cascade
c
c re-calculate wetted fractions, changed by cascade
c
      call fwetcal
c
c step vegetation canopy temperatures implicitly
c and calculate sensible heat and moisture fluxes
c
      call canopy
c
c step intercepted h2o due to evaporation
c
      call cascad2
c
c arbitrarily set veg temps & intercepted h2o for no-veg locations
c
      call noveg
c
c set net surface heat fluxes for soil and snow models
c
      do 110 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
        heatg(i) = solg(i) + firg(i) - fseng(i) -
     >             hvasug(i)*fvapg(i)
c
        heati(i) = soli(i) + firi(i) - fseni(i) -
     >             hvasui(i)*fvapi(i)
c
 110  continue
c
c step snow model
c
      call snow
c
c step soil model
c
      call soilctl
c
c return to main program
c
      return
      end
c 
