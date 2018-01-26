c Author: Bill Sacks
c Creation date: 02.10.10

c This module contains higher-level wrappers to the routines in
c ies-io.f. It was created because of a need for some features that
c require an explicit interface - e.g., optional subroutine arguments.
c 
c Note that this file contains some syntax / features that are not
c supported by fortran 77.

c ---------------------------------------------------------------------
      module netcdf_wrappers
c ---------------------------------------------------------------------
c 
      implicit none
      save

      contains

c ---------------------------------------------------------------------
      subroutine write4dvar (dat, nlev, mstep, filen, varname, ftime,
     >  tweight, tdate, scaling)
c ---------------------------------------------------------------------
c 
c writes a single time slice of a 4-d variable (lon x lat x level x
c time) to a netcdf file
c 
c inifile and inivar should already have been called for this netcdf
c file and output variable
c 
      use comgrid
      use compar
      use comwork
c
c Arguments
c 
      integer nlev              ! number of levels in dat (all levels will be written)

c the data array to write; it will be converted to a nlonsub x nlatsub x
c nlev array before being written
      real dat(npoi, nlev)  

      integer mstep             ! time slice index to write

      character*(*) 
     >  filen,                  ! name of output file (must already have been initialized with inifile)
     >  varname,                ! output variable name (must already have been initialized with inivar)
     >  tdate                   ! character date for time step

      real
     >  ftime(*),               ! real form of number of days run since iyear0
     >  tweight(*)              ! number of days in the time step

c scaling factor: all data items are multiplied by this factor
c (this avoids the creation of the possibly large temporary variable
c that would be created by doing the scaling inline in the function
c call - e.g., write4dvar(data/100., ...))
      real, optional :: scaling
c 
c Local variables
c 
      integer i, j, istat
      integer istart(4), icount(4) ! for writing vars

c if 'scaling' is provided, this array will hold one slice of the scaled
c array
c Note: we 'save' this variable so that it's not allocated on the stack
      real, save :: dat_scaled(npoi)
c 
      data istart / 1,1,1,1 /,
     >  icount / nlon,nlat,1,1 /
      istart(4) = mstep
      icount(1) = nlonsub
      icount(2) = nlatsub
      icount(3) = nlev

      if (nlev .gt. max3d) then
        write(*,*) 'ERROR IN WRITE4DVAR:'
        write(*,*)
     >    'variable has more levels than is accommodated by cdummy'
        write(*,*) 'nlev: ', nlev
        write(*,*) 'max3d: ', max3d
        stop 1
      end if

c copy data into cdummy, translating from a spatial dimension of npoi to
c a nlonsub x nlatsub grid
      do i = 1, nlev
        if (present(scaling)) then
          do j = 1, npoi
            dat_scaled(j) = dat(j,i) * scaling
            call vec2arr (dat_scaled, cdummy((i-1)*nlonsub*nlatsub + 1))
          end do
        else                    ! no scaling factor            
          call vec2arr (dat(1,i), cdummy((i-1)*nlonsub*nlatsub + 1))
        end if
      end do

c now write cdummy to the netcdf file
      call writevar(filen, varname, istart, icount, cdummy, ftime,
     >  tweight, tdate, istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in write4dvar, ', trim(varname)
        stop 1
      end if

      return
      end subroutine write4dvar

! --write4dvar_ver2: add support for user-specified istart and icount
!   write the chunk date <dat> all at once
!
! dimension format (lon,lat,level,time)
! dat: data array
! vardim: dimension of the variable (npoi,nlev,time), if only two dims,
!         let nlev = 1
c ---------------------------------------------------------------------
      subroutine write4dvar_ver2(filen,varname,vardim,istart,icount,dat,ftime,tweight,tdate,istat)
c ---------------------------------------------------------------------
c 
c writes all time slice of a 4-d variable (lon x lat x level x time) to a netcdf file
c 
c inifile and inivar should already have been called for this netcdf file and output variable
c 
      use comgrid
      use compar
      use comwork

c the data array to write; it will be converted to a nlonsub x nlatsub x nlev x ntime array before being written
      integer:: ntots,vardim(3)
      real   :: dat(vardim(1), vardim(2), vardim(3)) !usually vardim(1) = npoi
      real,dimension(:),allocatable:: varArray

      character*(*) 
     >  filen,                  ! name of output file (must already have been initialized with inifile)
     >  varname,                ! output variable name (must already have been initialized with inivar)
     >  tdate                   ! character date for time step

      real
     >  ftime(*),               ! real form of number of days run since iyear0
     >  tweight(*)              ! number of days in the time step

      integer i, j,idx, istat
      integer istart(4), icount(4) ! for writing vars
       
      !ntots = vardim(1)*vardim(2)*vardim(3)
      ntots = nlonsub*nlatsub*vardim(2)*vardim(3)  !use nlonsub*nlatsub instead of npoi to accommodate grid using surta, which could results in npoi<nonsub*nlatsub
      allocate(varArray(ntots),STAT=istat)

c copy data into varArray, translating from a spatial dimension of npoi to a nlonsub x nlatsub grid
      idx = 0
      do j = 1, vardim(3)       !time
        do i = 1, vardim(2)     !lev 

          idx = idx + 1
          call vec2arr (dat(:,i,j), varArray((idx-1)*nlonsub*nlatsub + 1) )

        end do !on i
      end do !on j

c now write varArray to the netcdf file
      call writevar(filen, varname, istart, icount, varArray, ftime, tweight, tdate, istat)
      if (istat .ne. 0) then
        write(*,*) 'ERROR in write4dvar_ver2, ', trim(varname)
        stop 1
      end if

      deallocate(varArray,STAT=istat)

      return
      end subroutine write4dvar_ver2

      end module netcdf_wrappers


      
