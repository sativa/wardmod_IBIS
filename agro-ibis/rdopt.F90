!-*-F90-*-
! ----------------------------------------------------------------------------! Module rdopt:  Contains an array of opened netCDF units/ncids, and routines
!                to log them, close them, and look for ncids that match a
!                file's pathname.
! rdopt_mmap_bufrsize is available for use with NF__OPEN using cmode=NF_MMAP.
!
! S. Slavin, 27 May 2016
! ----------------------------------------------------------------------------

module rdopt
  implicit none
  public

! Set to the default Lustre blocksize of 4 MB. 
  integer :: rdopt_mmap_bufrsize = 4*1024*1024  ! NF_MMAP buffer size in bytes

! Size of the ncid_map array/table of logged ncids that are open.
  integer, private, parameter :: max_units = 1000

  integer, private :: ncid_map( max_units ) = 0
  integer, private :: num_open = 0

contains

! ----------------------------------------------------------------------------

  subroutine log_open_ncid( ncid )
    implicit none

    integer, intent(in) :: ncid

    num_open = num_open + 1
    if ( num_open > max_units ) then
       stop '[rdopt.F90](log_open_ncid) Memory for ncid_map exceeded. ' &
            //'Increase max_units.'
    else
       ncid_map( num_open ) = ncid
    end if

  end subroutine log_open_ncid

! ----------------------------------------------------------------------------

  subroutine ncid_from_fname( fname, ncid )
    implicit none

    include 'netcdf.inc'

    character(len=*), intent(in) :: fname
    integer, intent(inout) :: ncid

    integer :: path_len
    character(len=NF_MAX_NAME+1) :: path_name

    integer :: i = 0
    integer :: j = 0
    integer :: jerr = 0

    path_name = repeat( ' ', len(path_name) )
    path_len = 0

    ncid = 0

! Loop through the files we've logged as being open
    do i = 1, num_open

! Grab the ncid for the i_th open file we've logged.
       j = ncid_map(i)

! Ask netCDF what the pathname for the open file is.
       jerr = nf_inq_path( j, path_len, path_name )
       if ( jerr /= NF_NOERR ) then
          stop '[rdopt.F90](ncid_from_fname): nf_inq_path() reported an error.'
       end if

! Check to see if it matches the filename passed into this routine.
       if ( trim(path_name) == trim(fname) ) then
          ncid = j
          exit
       end if
    end do
! Return with the ncid value of the open file matching the fname variable.
  end subroutine ncid_from_fname

! ----------------------------------------------------------------------------

  subroutine close_logged_ncids
    implicit none

    include 'netcdf.inc'

    integer :: jerr
    integer :: i

    do i = 1, num_open

       jerr = NF_CLOSE( ncid_map(i) )

       if ( jerr /= NF_NOERR ) then   ! Say something about error codes
          write(*,'(a)') '[rdopt.F90](close_logged_ncids) ' &
               // 'NF_CLOSE returned an error: ' &
               // NF_STRERROR(jerr)
       end if
    end do

    num_open = 0   ! Reset the open file counter
    ncid_map = 0   ! Reinitialize the ncid map array

  end subroutine close_logged_ncids

! ----------------------------------------------------------------------------

end module rdopt
