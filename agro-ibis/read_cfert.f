c Bill Sacks, Jan. 2010

c----------------------------------------------------------
      subroutine read_cfert
c----------------------------------------------------------
c 
c This subroutine reads in the text file that is used to set the cfert*
c variables. Until Jan. 2010, this code was in crops.f : planting.
c
      use comgrid
      use compar
      use comcrop
c
      implicit none
c
c Local variables
      integer 
     >  nread,       ! number of years read from file
     >  iy           ! current year being read
 
      character*100 fertdata,
     >              header

      real 
     >  this_cfertmaize,
     >  this_cfertsoy

      parameter (fertdata = 'hist.fert.1945.1996')  ! data file that contains historical
                                                    ! changes in fertilizer usage across upper Midwest US  

      cfert_year1 = 9999
      cfert_lastyear = 0
      nread = 0
      open(15, file = fertdata) 
      read(15, *) header

c Read until we reach end of file
      do
        read(15, *, end=20) iy, this_cfertmaize,
     >    this_cfertsoy
        if (iy .gt. max_cfert_years) then
          write(*,*) 'ERROR IN READ_CFERT'
          write(*,*) 'while reading ', trim(fertdata)
          write(*,*) 'year ', iy
          write(*,*) 'greater than max_cfert_years: ', max_cfert_years
          stop 1
        end if

        nread = nread + 1
        if (iy .lt. cfert_year1) then
          cfert_year1 = iy
        end if
        if (iy .gt. cfert_lastyear) then
          cfert_lastyear = iy
        end if

c file has quantities in kg/ha, so these are changed to kg/m2
        cfertmaize(iy) = this_cfertmaize * 1e-04
        cfertsoy(iy)   = this_cfertsoy   * 1e-04
c 
c temporarily assign wheat fertilizer usage to be 3 times that of
c soybeans
        cfertwheat(iy) = cfertsoy(iy) * 3.0
      enddo
 20   continue 
      close(15)

c Make sure that we read at least one year of data
      if (nread .eq. 0) then
        write(*,*) 'ERROR IN READ_CFERT'
        write(*,*) 'while reading ', trim(fertdata)
        write(*,*) 'No data read in'
        stop 1
      end if

c Make sure that we read data for every year between cfert_year1 and
c cfert_lastyear
      if (nread .ne. (cfert_lastyear - cfert_year1 + 1)) then
        write(*,*) 'ERROR IN READ_CFERT'
        write(*,*) 'while reading ', trim(fertdata)
        write(*,*) 'Some missing years'
        stop 1
      end if

      return
      end
      

