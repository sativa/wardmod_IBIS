c This file contains utilities for handling climate files
c Created by Bill Sacks, Jan. 2010
c
c Note that this file contains some syntax / features that are not
c supported by fortran 77.

c ---------------------------------------------------------------------
      module clim_file_utils
c ---------------------------------------------------------------------

      implicit none
      save

      contains

c ---------------------------------------------------------------------
c overlap
c
c     Given the start years of monthly and daily climate files, and the
c     number of years in monthly and daily files, return the start and
c     end years of overlap.
c
c     If nanom = 0, just return start and end years in daily files; if
c     nanomd = 0, just return start and end years in monthly files. If
c     both nanom and nanomd = 0, results are meaningless.
c ---------------------------------------------------------------------

      subroutine overlap(overlap_start, overlap_end, istyrm, nanom,
     >     istyrd, nanomd)

c Arguments

      integer 
     >     overlap_start,  ! output: start year of overlap
     >     overlap_end,    ! output: end year of overlap
     >     istyrm,         ! start year in monthly files
     >     nanom,          ! # years in monthly files
     >     istyrd,         ! start year in daily files
     >     nanomd          ! # years in daily files

c Local variables

      integer endyrm, endyrd

      endyrm = istyrm + nanom - 1
      endyrd = istyrd + nanomd - 1

      if (nanom .eq. 0) then
        overlap_start = istyrd
        overlap_end = endyrd
      else if (nanomd .eq. 0) then
        overlap_start = istyrm
        overlap_end = endyrm
      else
        overlap_start = max(istyrm, istyrd)
        overlap_end = min(endyrm, endyrd)
      end if

      return
      end subroutine overlap

c ---------------------------------------------------------------------
c year_to_read
c
c     Determine the year index to read in a climate file. This function
c     handles cycling through years, and ensures that the monthly and
c     daily files stay in sync by only using overlapping years in the
c     cycling. (Note that, if there are no daily files, then
c     overlap_start and overlap_end will just give the start and end of
c     the monthly file, and similarly if there are no monthly files.)
c
c     Returns an integer giving the index of the year to read, where 0 =
c     1st year, 1 = 2nd year, etc.
c
c     If iyear < istyr, a negative value is returned; it is up to the
c     calling routine to check this error condition.
c ---------------------------------------------------------------------

      integer function year_to_read(iyear, istyr, overlap_start,
     >     overlap_end)
      
c Arguments
      
      integer 
     >     iyear,          ! calendar year desired
     >     istyr,          ! start year in this file
     >     overlap_start,  ! start year of overlap between monthly & daily files
     >     overlap_end     ! end year of overlap between monthly & daily files

c Local variables

      integer noverlap, years_past_end

      if (iyear .le. overlap_end) then
        year_to_read = iyear - istyr
      else
        noverlap = overlap_end - overlap_start + 1
        years_past_end = iyear - overlap_end
        year_to_read = (overlap_start - istyr) +
     >       mod((years_past_end - 1), noverlap)
      end if

      return
      end function year_to_read


c ---------------------------------------------------------------------
c climate_full_fname
c
c     Given a file name (fname) without the .nc extension, and without
c     the possible _YYYY suffix, return the full path and file name of
c     this input climate file in the variable full_fname (full_fname
c     must be large enough to hold the result).
c
c Examples (assuming base_path = 'input/')
c  fname = 'xyz', year = 1901, use_year = true  -> full_fname = 'xyz_1901.nc'
c  fname = 'xyz', year = 123, use_year = true   -> full_fname = 'xyz_0123.nc'
c  fname = 'xyz', year = 0, use_year = true     -> full_fname = 'xyz_0000.nc'
c  fname = 'xyz', year = 1901, use_year = false -> full_fname = 'xyz.nc'
c ---------------------------------------------------------------------

      subroutine climate_full_fname(full_fname, fname, year, use_year)

c Arguments
 
c     full_fname: Will hold the result
      character*(*) full_fname

c     fname: file name without the .nc extension, and without the
c     possible _YYYY suffix
      character*(*) fname

c     year: specifies the _YYYY suffix, if use_year = true
      integer year

c     use_year: If true, there is a _YYYY suffix; if false, no _YYYY
c     suffix (see examples in the header)
      logical use_year

c Local variables
c     base_path: constant giving the path to all the climate input files
      character*64 base_path
c     suffix: holds optional _YYYY, plus .nc
      character*8 suffix

      parameter(base_path = 'input/')

c External
      integer lenchr


      if (use_year) then
        if (year .lt. 0 .or. year .gt. 9999) then
          write(*,*) 'ERROR in climate_full_fname:'
          write(*,*) 'year must be >= 0 and <= 9999'
          stop 1
        else
c         Set suffix = '_YYYY.nc'
          write(suffix, '("_", i4.4, ".nc")') year 
        end if
      else  ! .not. use_year
        suffix = ".nc"
      end if

      if (lenchr(base_path) + lenchr(fname) + lenchr(suffix) .gt.
     >     len(full_fname)) then
        write(*,*) 'ERROR in climate_full_fname:'
        write(*,*) 'full_fname too small to hold result'
        stop 1
      end if

      full_fname = trim(base_path) // trim(fname) // trim(suffix)
      return
      end subroutine climate_full_fname

c ---------------------------------------------------------------------
c climate_years_avail
c
c     Given the start and end years of available climate data, allocate
c     and return two vectors: one containing a list of available leap
c     years, and the other containing a list of available non-leap
c     years.
c ---------------------------------------------------------------------

      subroutine climate_years_avail(avail_leaps, avail_nonleaps,
     >  first_year, last_year)

c Arguments

c These arrays will hold lists of years of available climate data,
c divided into leap years and non-leap years. Available years are
c defined as all years between first_year and last_year, inclusive. If
c the arrays are already allocated, they will be deallocated. They will
c then be allocated to be exactly the necessary size to hold their
c respective lists (or size 0 if there are no years of that type).
      integer, dimension(:), allocatable, intent(inout) :: avail_leaps
      integer, dimension(:), allocatable, intent(inout) ::
     >  avail_nonleaps

c First and last year of available climate data.
      integer, intent(in) :: first_year
      integer, intent(in) :: last_year

c Local variables

      integer n_leap, n_nonleap ! number of leap years and non-leap years
      integer yr
      integer i_l, i_nl  ! indices into avail_leaps and avail_nonleaps

c Externals

      logical is_leap        ! Function: is the given year a leap year?


      if (allocated(avail_leaps)) then
        deallocate(avail_leaps)
      end if
      if (allocated(avail_nonleaps)) then
        deallocate(avail_nonleaps)
      end if

c First pass: determine number of leap years and non-leap years
      n_leap = 0
      n_nonleap = 0
      do yr = first_year, last_year
        if (is_leap(yr)) then
          n_leap = n_leap + 1
        else
          n_nonleap = n_nonleap + 1
        end if
      end do

      allocate(avail_leaps(n_leap))
      allocate(avail_nonleaps(n_nonleap))

c Second pass: fill the vectors with the appropriate years
      i_l = 1
      i_nl = 1
      do yr = first_year, last_year
        if (is_leap(yr)) then
          avail_leaps(i_l) = yr
          i_l = i_l + 1
        else
          avail_nonleaps(i_nl) = yr
          i_nl = i_nl + 1
        end if
      end do

      end subroutine climate_years_avail
      end module clim_file_utils
      
        

c Testing code: uncomment and compile this file to test some of the
c subroutines contained here. If all goes well, you should not see any
c instances of 'ERROR'.
c
c Compilation example:
c     ifort -g -O0 -warn all -traceback -fpe0 -check all -132 clim_file_utils.f utilities.f math.f
c
c$$$      program testing
c$$$
c$$$      use clim_file_utils
c$$$      integer overlap_start, overlap_end
c$$$      integer, dimension(:), allocatable :: avail_leaps, avail_nonleaps
c$$$      
c$$$      if (year_to_read(1902, 1901, 1948, 2000) .ne. 1) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      if (year_to_read(2000, 1901, 1948, 2000) .ne. 99) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      if (year_to_read(2001, 1901, 1948, 2000) .ne. 47) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      if (year_to_read(2055, 1901, 1948, 2000) .ne. 48) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      if (year_to_read(1948, 1948, 1948, 2000) .ne. 0) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      if (year_to_read(2055, 1948, 1948, 2000) .ne. 1) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      if (year_to_read(1955, 1901, 1901, 2000) .ne. 54) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      if (year_to_read(2055, 1901, 1901, 2000) .ne. 54) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      if (year_to_read(2155, 1901, 1901, 2000) .ne. 54) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$c Check error in year_to_read (trying to read before start of file)
c$$$c (Note: we don't care what the exact value is - we just want to make
c$$$c sure that the return value is negative)
c$$$      if (year_to_read(1900, 1901, 1948, 2000) .ge. 0) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      call overlap(overlap_start, overlap_end,1901,100,1948,55)
c$$$      if (overlap_start .ne. 1948 .or. overlap_end .ne. 2000) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      call overlap(overlap_start, overlap_end,1948,55,1901,100)
c$$$      if (overlap_start .ne. 1948 .or. overlap_end .ne. 2000) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      call overlap(overlap_start, overlap_end, 1901,100,1948,0)
c$$$      if (overlap_start .ne. 1901 .or. overlap_end .ne. 2000) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      call overlap(overlap_start, overlap_end, 1901,0,1948,55)
c$$$      if (overlap_start .ne. 1948 .or. overlap_end .ne. 2002) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$      call climate_years_avail(avail_leaps, avail_nonleaps, 1901, 1910)
c$$$      if (size(avail_leaps) .ne. 2 .or. size(avail_nonleaps) .ne. 8
c$$$     >  .or. .not. all(avail_leaps .eq. (/1904, 1908/))
c$$$     >  .or. .not. all (avail_nonleaps .eq.
c$$$     >  (/1901, 1902, 1903, 1905, 1906, 1907, 1909, 1910/))) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$
c$$$c Check climate_years_avail with one array being size 0
c$$$      call climate_years_avail(avail_leaps, avail_nonleaps, 1901, 1903)
c$$$      if (size(avail_leaps) .ne. 0 .or. size(avail_nonleaps) .ne. 3
c$$$     >  .or. .not. all(avail_nonleaps .eq. (/1901, 1902, 1903/))) then
c$$$        print *, 'ERROR'
c$$$      end if
c$$$      
c$$$      stop0
c$$$      end
