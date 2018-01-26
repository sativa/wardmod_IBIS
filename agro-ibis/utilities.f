c
c #    #   #####     #    #          #     #####     #    ######   ####
c #    #     #       #    #          #       #       #    #       #
c #    #     #       #    #          #       #       #    #####    ####
c #    #     #       #    #          #       #       #    #            #
c #    #     #       #    #          #       #       #    #       #    #
c  ####      #       #    ######     #       #       #    ######   ####
c
c ---------------------------------------------------------------------
      subroutine scopy (nt, arr, brr)
c ---------------------------------------------------------------------
c
c copies array arr to brr,for 1st nt words of arr
c
      implicit none
c
c Arguments
c
      integer nt     
      real arr(nt),    ! input
     >     brr(nt)     ! output
c
c Local variables
c
      integer ia
c
      do 100 ia = 1, nt
        brr(ia) = arr(ia)
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine const (arr, nar, value)
c ---------------------------------------------------------------------
c
c sets all elements of real vector arr to value
c
      implicit none
c
c Arguments
c
      integer nar
c     
      real value
      real arr(nar)
c
c Local variables
c
      integer j
c
      do 100 j = 1, nar
        arr(j) = value
 100  continue
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine const_int (arr, nar, value)
c ---------------------------------------------------------------------
c
c sets all elements of int vector arr to value
c
      implicit none
c
c Arguments
c
      integer nar
c     
      integer value
      integer arr(nar)
c
c Local variables
c
      integer j
c
      do 100 j = 1, nar
        arr(j) = value
 100  continue
c
      return
      end
c
c
c
c Set phenology onset flags
c intially to .FALSE. at beginning of first year
c and phenology offset flags
c intially to .TRUE. at beginning of first year
c
c ---------------------------------------------------------------------
      subroutine logict (arrl, nar)
c ---------------------------------------------------------------------
      logical arrl(nar)
c
      do 100 j = 1, nar
        arrl(j) = .TRUE.
 100  continue
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine logicf (arrl, nar)
c ---------------------------------------------------------------------
      logical arrl(nar)
c
      do 100 j = 1, nar
        arrl(j) = .FALSE.
 100  continue
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine endrun
c ---------------------------------------------------------------------
c
c stops gracefully
c
      stop0
      end
c
c 
c ---------------------------------------------------------------------
      real function cvmgt (x,y,l)
c ---------------------------------------------------------------------
c
c chooses between two things.  Used in canopy.f
c
      implicit none
c
      logical l
      real x, y
c
      if (l) then
        cvmgt = x
      else
        cvmgt = y
      endif
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine arr2vec(array, vect)
c ---------------------------------------------------------------------
c
c extracts land points from array and puts them into a land-only vector
c
      use comgrid
      use compar
      use combcs
      use comcrop
      use comsoi
      use comnitr
c
      implicit none
c
c Arguments
c
      real array(nlonsub,nlatsub), vect(npoi)
c
c Local variables
c
      integer j,i,npts
c
      npts = 0
c
      do 10 j = 1, nlatsub
        do 20 i = 1, nlonsub
          if (lmask(i,j) .eq. 1) then
            npts = npts + 1
            vect(npts) = array(i,j)
          end if
 20     continue
 10   continue
c
      if (npts .ne. npoi) then
        write (*,*) 'ERROR in arr2vec'
        write (*,*) 'npts not equal to npoi'
        write (*,*) 'npts = ', npts, ' npoi = ', npoi
        stop 1
      end if
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vec2arr(vect, array)
c ---------------------------------------------------------------------
c
c puts vector of land-only points back into array
c
      use comgrid
      use compar
      use comwork
c
      implicit none
c
c Arguments
c
      real array(nlonsub,nlatsub), vect(npoi)
c
c Local variables
c
      integer j, i, ii, jj
c
      do 10 j = 1, nlatsub
        do 20 i = 1, nlonsub
          array(i,j) = OCEAN
 20     continue
 10   continue
c
      do 30 i = 1, npoi
        ii = lonindex(i)
        jj = latindex(i)
        array(ii,jj) = vect(i)
 30   continue
c
      return
      end
c
c ---------------------------------------------------------------------
      subroutine vec2arr_ver2(vect, array)
c ---------------------------------------------------------------------
c
c puts all points (land and water to array)
c land variable will be from vect, water set as 0 for now
c
      use comgrid
      use compar
      use comwork
c
      implicit none
c
c Arguments
c
      real array(nlonsub,nlatsub), vect(npoi)
c
c Local variables
c
      integer j, i, ii, jj
c
      array = 0.0 

c
      do 30 i = 1, npoi
        ii = lonindex(i)
        jj = latindex(i)
        array(ii,jj) = vect(i)
 30   continue
c
      return
      end subroutine vec2arr_ver2

c
c
c ---------------------------------------------------------------------
c lenchr - find index of last non-blank, non-null
c ---------------------------------------------------------------------
c
      integer function lenchr (ch)
c
c returns position of last non-blank,null character in ch,
c or 1 if ch is all blanks
c
      implicit none
c
c Arguments
c
      character*(*) ch
c
      integer i
c
      do i = len(ch), 1, -1
        if (ch(i:i).ne.' '.and.ch(i:i).ne.char(0)) then
           lenchr = i
           return
        endif
      enddo
c
      lenchr = 1
c
      return
      end
c
c ---------------------------------------------------------------------
      real function get_time(istep)
c ---------------------------------------------------------------------
c
c Return the current time of day since midnight (in seconds)
c
      use compar  ! needed for dtime
c
      implicit none
c
c Arguments
c
      integer istep  ! timestep counter (per day)

      get_time = (istep - 1) * dtime

      return
      end
c
c ---------------------------------------------------------------------
      logical function is_leap(yr)
c ---------------------------------------------------------------------
c
c Return true if yr is a leap year, false otherwise
c
      implicit none
c
c Arguments
c
      integer yr

c leap years occur every calendar year that is divisible by 4 (such as
c 1980, 1984, etc.) *except* for years which begin centuries not divisible
c by 400
c
c for example, 1700, 1800, and 1900 are *not* leap years, but 1600
c and 2000 *are* leap years
      
      if (mod(yr,4) .eq. 0) then
        if (mod(yr,100) .ne. 0) then
          is_leap = .TRUE.
        else if (mod(yr/100,4) .eq. 0) then
          is_leap = .TRUE.
        else  ! year is divisible by 100, but not by 400
          is_leap = .FALSE.
        end if
      else  ! year not divisible by 4
        is_leap = .FALSE.
      end if

      return
      end
c
c ---------------------------------------------------------------------
      subroutine random_choice(val, arr, n, seed)
c ---------------------------------------------------------------------
c
c Given an array of integers (arr), randomly choose and return an
c element between arr(1) and arr(n), inclusive. The random element is
c returned in 'val'. 
c
      implicit none
c
c Arguments
c 
      integer val      ! output: the randomly-chosen value
      integer arr(*)   ! input: array of values to choose from
      integer n        ! input: we just consider the first n elements of arr (n must be > 0)
      integer seed     ! input / output: random number seed (modified in call to random number generator)
c
c Local variables
c
      integer i
c
c Externals
c 
      real ran2        ! Function: random number generator


      if (n .le. 0) then
        write(*,*) 'ERROR in random_choice: n must be > 0'
        stop 1
      end if

      i = int(ran2(seed)*n) + 1 ! random number between 1 and n, inclusive
      val = arr(i)

      return
      end

c ---------------------------------------------------------------------
      integer function random_seed()
c ---------------------------------------------------------------------
c 
c Return a negative integer that can be used to seed the random number
c generator. This seed is based on the current CPU time. The algorithm
c used here could probably be improved (e.g., using seeds that differ
c more widely for two nearby times, and spanning the range of negative
c integers more fully), but it should be sufficient for our purposes.
c
      implicit none
c
c Local variables
c
      integer 
     >  seconds,         ! # seconds since start of year (approx.)
     >  seed             ! the seed to return

      integer values(8)  ! holds values from date_and_time

c Arbitrary positive integer for scaling the year, so that two runs done
c one year apart are less likely to end up with the same seed. The exact
c value is not important, but if you change this, it should be greater
c than 0, and less than about 1 million to avoid integer overflow.
      integer year_scale
      parameter (year_scale = 7919)  ! 7919 happens to be the 1000th prime


      call date_and_time(values=values)

c Determine approximate seconds since start of year (but assuming that
c all months have 31 days for simplicity)
      seconds = values(7)    ! seconds
     >  + 60*(values(6)      ! minutes
     >  + 60*(values(5)      ! hours
     >  + 24*((values(3)-1)  ! day of month
     >  + 31*((values(2)-1)  ! month
     >  ))))      

c Create the seed by adding the following to the seconds:
c - current milliseconds, so that two runs done in quick succession have
c   different seeds
c - current year, scaled by some arbitrary positive integer, so that two
c   runs done exactly one year apart are less likely to have the same
c   seed
c Then multiplying by negative one to make it negative

      seed = seconds
     >  + values(8)             ! milliseconds
     >  + year_scale*values(1)  ! year

      seed = -1*seed

c At this point, seed should generally be negative, unless there was
c integer overflow. Here we ensure that the final value of seed truly is
c negative:
      if (seed .gt. 0) then
        seed = -1*seed
      else if (seed .eq. 0) then
        seed = -1
      end if

      random_seed = seed
      return
      end

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : find possible neighbor cells index in global 1D
! 
! 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine find_neighbors(idx,n_neighbors,idx_neighbors)

      use comgrid
      use comwork

      implicit NONE

      !input
      integer:: idx  !global 1D index (1~npoi)
      !output
      integer:: n_neighbors,idx_neighbors(8)  !number of neighbors, and their global index
      !local
      integer i,ii,jj,idx_ii,idx_jj
      logical flg_neighbor

      n_neighbors   = 0
      idx_neighbors = 0

      idx_ii = lonindex(idx)
      idx_jj = latindex(idx)
      do i = 1, npoi
        flg_neighbor = .false.
        ii = lonindex(i)
        jj = latindex(i)

        !counter-clock wise searching
        if     ( ii == idx_ii   .and. jj == idx_jj-1) then 
          flg_neighbor = .true.
        elseif ( ii == idx_ii+1 .and. jj == idx_jj-1) then
          flg_neighbor = .true.
        elseif ( ii == idx_ii+1 .and. jj == idx_jj)   then
          flg_neighbor = .true.
        elseif ( ii == idx_ii+1 .and. jj == idx_jj+1) then
          flg_neighbor = .true.
        elseif ( ii == idx_ii   .and. jj == idx_jj+1) then
          flg_neighbor = .true.
        elseif ( ii == idx_ii-1 .and. jj == idx_jj+1) then
          flg_neighbor = .true.
        elseif ( ii == idx_ii-1 .and. jj == idx_jj)   then
          flg_neighbor = .true.
        elseif ( ii == idx_ii-1 .and. jj == idx_jj-1) then
          flg_neighbor = .true.
        end if

        if(flg_neighbor) then
          n_neighbors = n_neighbors + 1
          idx_neighbors(n_neighbors) = i
        end if


      end do !on i = 1,npoi
      

      end subroutine  find_neighbors
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! purpose : error handler for commandline executioin
! 
! 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine error_handler_command(exit_s,cmd_s,message, note)
        implicit NONE

        integer:: exit_s, cmd_s
        character*200 message, note

        if( cmd_s /= 0 ) then
          write(6,*) trim(note)
          write(6,*) '(error_handler_command): '//trim(message)
          if(cmd_s == -1) then
            write(6,'(A)') '(error_handler_command): command execution not supported on this system'
            write(6,'(A)') 'Stopping the program...'
            stop
          end if
        end if

        if(exit_s /= 0 ) then
          write(6,'(A)') '(error_handler_command): command execution failed'
          write(6,'(A)') 'Stopping the program...'
          stop
        end if

      
        return
      end subroutine
