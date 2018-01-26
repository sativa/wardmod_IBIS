c This file contains functions for calculating GDD in a few different
c ways
c 
c Created by Bill Sacks, 03.24.10

c ---------------------------------------------------------------------
      real function gdd_natveg(tmean, base)
c ---------------------------------------------------------------------
c
c Calculate one day's GDD using the function that applies to natural
c vegetation
c
      implicit none
c
c Arguments
c
c Note: units don't matter, as long as the units of tmean and base are
c the same
      real tmean                ! mean temperature on the given day
      real base                 ! base temperature for GDD calculations

      gdd_natveg = max(0.0, (tmean - base))
      
      return
      end


c ---------------------------------------------------------------------
      real function gdd_soil(tsoi, base, gdd_max)
c ---------------------------------------------------------------------
c
c Calculate one day's GDD using the function that applies to soil
c temperatures
c
c WJS (03.24.10): I am copying the meaning of gdd_max that is used in
c the code, although I think that it is often applied in a slightly
c different way: if tmean > gdd_max, then tmean is set equal to gdd_max,
c then base is subtracted from this. e.g., if tmean = 35, base = 10,
c gdd_max = 30, the current function will give gdd = 25, but according
c to other sources, gdd would be capped at 20 with these values of base
c and gdd_max.
c
      implicit none
c
c Arguments
c
c Note: units don't matter, as long as the units of tsoi, base and
c gdd_max are the same
      real tsoi                 ! mean soil temperature on the given day
      real base                 ! base temperature for GDD calculations
      real gdd_max              ! maximum GDD that can be accumulated
                                ! per day (must be > 0)

      gdd_soil = max(0.0, min((tsoi - base), gdd_max))
      
      return
      end



c ---------------------------------------------------------------------
      real function gdd_crop(tmean, tmax, tmin, base, gdd_max)
c ---------------------------------------------------------------------
c
c Calculate one day's GDD using the function that applies to crops. This
c is a wrapper for one of the specific gdd_crop_* functions. Some of the
c parameters may be unused depending on the specific gdd_crop_* function
c that is being used: we include all possible parameters so that the
c call to gdd_crop doesn't have to change if the GDD implementation
c changes.
c
c The point of this wrapper is that the GDD calls throughout the code
c can simply be to gdd_crop(tmean, tmax, tmin, base, gdd_max), and then
c the only thing we need to change in order to change the GDD
c implementation is the body of this function - we don't have to change
c any of the function calls in the rest of the code.
c
      implicit none
c
c Arguments
c
c Note: units don't matter, as long as the units of tmean, tmax, tmin,
c base and gdd_max are the same
      real tmean                ! mean temperature on the given day
      real tmax                 ! daily max temperature
      real tmin                 ! daily min temperature
      real base                 ! base temperature for GDD calculations
      real gdd_max              ! the meaning of this parameter differs
                                ! depending on the specific function
                                ! being used
c
c Externals
c
      real gdd_crop_tmean, gdd_crop_modified

c All but one of the following lines should be commented out:

c      gdd_crop = gdd_crop_tmean(tmean, base, gdd_max)

      gdd_crop = gdd_crop_modified(tmax, tmin, base, gdd_max)

      return 
      end

c ---------------------------------------------------------------------
c SPECIFIC FUNCTIONS CALLED BY GDD_CROP
c ---------------------------------------------------------------------

c NOTE: These specific GDD functions are not intended to be called
c directly. Instead, call the generic gdd_crop.

c ---------------------------------------------------------------------
      real function gdd_crop_tmean(tmean, base, gdd_max)
c ---------------------------------------------------------------------
c
c Calculate one day's GDD using one possible formulation over crops: use
c daily mean temperature, a base temperature of 'base', and allow no
c more than 'gdd_max' GDD to be accumulated in a day.
c
c NOTE: This is not intended to be called directly. Instead, call the
c generic gdd_crop.
c
c WJS (03.24.10): I am copying the meaning of gdd_max that is used in
c the code, although I think that it is often applied in a slightly
c different way: I think that many people calculate GDD as: if tmean >
c gdd_max, then tmean is set equal to gdd_max, and then base is
c subtracted from this. e.g., if tmean = 35, base = 10, gdd_max = 30,
c the current function will give gdd = 25, but according to other
c sources, GDD would be capped at 20 with these values of base and
c gdd_max.
c
      implicit none
c
c Arguments
c
c Note: units don't matter, as long as the units of tmean, base and
c gdd_max are the same
      real tmean                ! mean temperature on the given day
      real base                 ! base temperature for GDD calculations
      real gdd_max              ! maximum GDD that can be accumulated
                                ! per day (must be > 0)

      gdd_crop_tmean = max(0.0, min((tmean - base), gdd_max))
      
      return
      end

c ---------------------------------------------------------------------
      real function gdd_crop_modified(tmax, tmin, base, gdd_max)
c ---------------------------------------------------------------------
c
c Calculate one day's GDD using the "modified" GDD formula (e.g., see
c http: //www.agry.purdue.edu/ext/corn/news/timeless/HeatUnits.html)
c
c NOTE: This is not intended to be called directly. Instead, call the
c generic gdd_crop.
c
      implicit none
c
c Arguments
c
c Note: units don't matter, as long as the units of tmax, tmin, base and
c gdd_max are the same
      real tmax                 ! daily max temperature
      real tmin                 ! daily min temperature
      real base                 ! base temperature for GDD calculations
      real gdd_max              ! max temperature for GDD calculations
                                ! (must be greater than 'base') (note
                                ! the different meaning compared to
                                ! gdd_crop_tmean)
c
c Local variables
c
      real 
     >  tmax_adjusted,
     >  tmin_adjusted,
     >  tmean_adjusted

      tmax_adjusted = max(min(tmax, gdd_max), base)
      tmin_adjusted = max(min(tmin, gdd_max), base)

c WJS (03.30.10): Note that this calculation of tmean differs from the
c calculations used elsewhere in the model, which use the weighted mean
c (based on tmin_wt and tmax_wt). I don't think this inconsistency
c matters, but if it makes you feel better, then don't think about this
c quantity as the daily mean temperature, but instead think about it as
c just some empirical quantity based on tmax and tmin. After all, this
c GDD formula is just some empirical formula that helps us reproduce
c observations well in practice; tmean_adjusted could just as easily be
c defined as, say, tmax_adjusted*0.8 + tmin_adjusted*0.2, if tmax was
c found to affect the accumulation of heat units in the plant much more
c than tmin did.
c
c WJS (07.12.10): The above is no longer an issue, since tmin_wt and
c tmax_wt have been changed to 0.5.
c
      tmean_adjusted = (tmax_adjusted + tmin_adjusted) / 2.

      gdd_crop_modified = tmean_adjusted - base

      return
      end

c ---------------------------------------------------------------------
c END SPECIFIC FUNCTIONS CALLED BY GDD_CROP
c ---------------------------------------------------------------------


c ---------------------------------------------------------------------
c TESTING SUBROUTINES FOR THE ABOVE FUNCTIONS
c ---------------------------------------------------------------------

c Each subroutine will print 'ERROR' followed by a test ID if there are
c any errors detected.

c ---------------------------------------------------------------------
      subroutine test_gdd_natveg
c ---------------------------------------------------------------------
      implicit none

      type test
      real tmean                ! tmean argument
      real base                 ! base argument
      real expected             ! expected value
      character(len=64) id      ! identifier of the test
      end type test

      type(test), dimension(2), parameter :: tests = (/
     >  test(10.0, 3.0, 7.0, "gdd_natveg: simple"),
     >  test(10.0, 11.0, 0.0, "gdd_natveg: 0 gdd") /)
      
      real, parameter :: epsilon = 1.e-6
      integer i
      real gdd_natveg
      
      do i = 1, size(tests)
        if (abs(gdd_natveg(tests(i)%tmean, tests(i)%base) -
     >    tests(i)%expected) .gt. epsilon) then
          print *, 'ERROR: ', tests(i)%id
        end if
      end do

      return
      end 


c ---------------------------------------------------------------------
      subroutine test_gdd_soil
c ---------------------------------------------------------------------
      implicit none
      
      type test
      real tsoi                 ! tsoi argument
      real base                 ! base argument
      real gdd_max              ! gdd_max argument
      real expected             ! expected value
      character(len=64) id      ! identifier of the test
      end type test

      type(test), dimension(3), parameter :: tests = (/
     >  test(15.0, 10.0, 30.0, 5.0, "gdd_soil: simple"),
     >  test(9.0, 10.0, 30.0, 0.0, "gdd_soil: 0 gdd"),
     >  test(100.0, 10.0, 30.0, 30.0, "gdd_soil: max") /)

      real, parameter :: epsilon = 1.e-6
      integer i
      real gdd_soil

      do i = 1, size(tests)
        if (abs(gdd_soil(tests(i)%tsoi, tests(i)%base, tests(i)%gdd_max)
     >    -tests(i)%expected) .gt. epsilon) then
          print *, 'ERROR: ', tests(i)%id
        end if
      end do

      return
      end


c ---------------------------------------------------------------------
      subroutine test_gdd_crop_tmean
c ---------------------------------------------------------------------
      implicit none
      
      type test
      real tmean                ! tmean argument
      real base                 ! base argument
      real gdd_max              ! gdd_max argument
      real expected             ! expected value
      character(len=64) id      ! identifier of the test
      end type test

      type(test), dimension(3), parameter :: tests = (/
     >  test(15.0, 10.0, 30.0, 5.0, "gdd_crop_tmean: simple"),
     >  test(9.0, 10.0, 30.0, 0.0, "gdd_crop_tmean: 0 gdd"),
     >  test(100.0, 10.0, 30.0, 30.0, "gdd_crop_tmean: max") /)

      real, parameter :: epsilon = 1.e-6
      integer i
      real gdd_crop_tmean

      do i = 1, size(tests)
        if (abs(gdd_crop_tmean(tests(i)%tmean, tests(i)%base,
     >    tests(i)%gdd_max) -tests(i)%expected) .gt. epsilon) then
          print *, 'ERROR: ', tests(i)%id
        end if
      end do

      return
      end
      

c ---------------------------------------------------------------------
      subroutine test_gdd_crop_modified
c ---------------------------------------------------------------------
      implicit none
      
      type test
      real tmax                 ! tmax argument
      real tmin                 ! tmin argument
      real base                 ! base argument
      real gdd_max              ! gdd_max argument
      real expected             ! expected value
      character(len=64) id      ! identifier of the test
      end type test

      type(test), dimension(7), parameter :: tests = (/
     >  test(80., 55., 50., 86., 17.5, "gdd_crop_modified: simple"),
     >  test(90., 72., 50., 86., 29., "gdd_crop_modified: high max"),
     >  test(68., 41., 50., 86., 9., "gdd_crop_modified: low min"),
     >  test(90., 88., 50., 86., 36.,
     >  "gdd_crop_modified: high max and min"),
     >  test(43., 41., 50., 86., 0.,
     >  "gdd_crop_modified: low max and min"),
c Note for the following two tests: it shouldn't matter if max and min are
c reversed: we should still get the same results
     >  test(72., 90., 50., 86., 29., "gdd_crop_modified: high min"),
     >  test(41., 68., 50., 86., 9., "gdd_crop_modified: low max") /)



      real, parameter :: epsilon = 1.e-6
      integer i
      real gdd_crop_modified

      do i = 1, size(tests)
        if (abs(gdd_crop_modified(tests(i)%tmax, tests(i)%tmin,
     >    tests(i)%base, tests(i)%gdd_max) -tests(i)%expected) .gt.
     >    epsilon) then
          print *, 'ERROR: ', tests(i)%id
        end if
      end do

      return
      end


c ---------------------------------------------------------------------
      subroutine test_gdd_crop
c ---------------------------------------------------------------------
c
c Here we just do some simple tests of gdd_crop to make sure that the
c parameters are being passed correctly - more detailed tests of the
c possible underlying functions are done elsewhere.
c
      implicit none
      
      type test
      real tmean                ! tmean argument
      real tmax                 ! tmax argument
      real tmin                 ! tmin argument
      real base                 ! base argument
      real gdd_max              ! gdd_max argument
      real expected             ! expected value
      character(len=64) id      ! identifier of the test
      end type test

c Tests to run if gdd_crop uses gdd_crop_tmean
      type(test), dimension(1), target :: tests_tmean = (/
     >  test(15., 16., 11., 10., 30., 5., "gdd_crop: tmean") /)

c Tests to run if gdd_crop uses gdd_crop_modified
      type(test), dimension(1), target :: tests_modified = (/
     >  test(15., 16., 11., 10., 30., 3.5, "gdd_crop: modified") /)
      
c Pointer to the correct tests to run
      type(test), dimension(:), pointer :: tests

      real, parameter :: epsilon = 1.e-6
      integer i
      real gdd_crop

c Select the correct tests to run based on the implementation of
c gdd_crop:
      tests => tests_modified

      do i = 1, size(tests)
        if (abs(gdd_crop(tests(i)%tmean, tests(i)%tmax, tests(i)%tmin,
     >    tests(i)%base, tests(i)%gdd_max) - tests(i)%expected) .gt.
     >    epsilon) then
          print *, 'ERROR: ', tests(i)%id
        end if
      end do

      return
      end

      
c ---------------------------------------------------------------------
c END TESTING SUBROUTINES FOR THE ABOVE FUNCTIONS
c ---------------------------------------------------------------------


c Testing code: uncomment and compile this file to test the functions
c contained here. If all goes well, you should not see any instances of
c 'ERROR'.
c
c Compilation example:
c     ifort -g -O0 -warn all -traceback -fpe0 -check all -132 gdd.f
c
c$$$      program testing
c$$$      implicit none
c$$$
c$$$      call test_gdd_natveg
c$$$      call test_gdd_soil
c$$$      call test_gdd_crop_tmean
c$$$      call test_gdd_crop_modified
c$$$      call test_gdd_crop
c$$$
c$$$      stop
c$$$      end




