c This file contains functions similar to those contained in comsat.h -
c i.e., functions dealing with humidity conversions. Making them real
c functions (as opposed to statement functions) allows slightly more
c complex & easier-to-follow logic.
c
c Created by Bill Sacks, 02.12.10

c ---------------------------------------------------------------------
      real function tdew_from_e (e)
c ---------------------------------------------------------------------
c
c calculates dew point temperature (K) from vapor pressure (n/m**2);
c from Campbell & Norman 1998 ("Introduction to Environmental
c Biophysics"), eqn. 3.14.
c
      use compar
c
      implicit none
c
      include 'comsat.h'
c
c Arguments
c
      real e     ! vapor pressure (n/m**2)

c Use different coefficients depending on whether tdew is above freezing
c (in which case e >= tdew_a) or below freezing (in which case e <=
c tdew_a). (But note that the function is continuous around e = tdew_a.)
      if (e .ge. tdew_a) then
        tdew_from_e = (tdew_c_w * log(e/tdew_a)) /
     >   (tdew_b_w - log(e/tdew_a)) + 273.16
      else
        tdew_from_e = (tdew_c_i * log(e/tdew_a)) /
     >   (tdew_b_i - log(e/tdew_a)) + 273.16
      end if

      return
      end

c ---------------------------------------------------------------------
      real function e_from_tdew (tdew)
c ---------------------------------------------------------------------
c
c calculates vapor pressure (n/m**2) from dew point temperature (K);
c derived by inverting tdew_from_e
c
c not coincidentally, this is identical to the function esat(t) given by
c eqn. 3.8 in Campbell & Norman 1998
c
c WJS 02.12.10: In theory, this function should give the same result as
c esat(tdew). But, since they use different approximations, they differ
c slightly in practice. From some tests, the differences are < 0.2% for
c temperatures between 0-50 C, but the differences can be greater for
c temperatures less than 0 C: e.g., ~ 1.6% difference for tdew = -40 C.
c However, having consistency between e_from_tdew and tdew_from_e seems
c more important than having consistency between e_from_tdew and esat.
c At least as of now (02.12.10), the inconsistencies between e_from_tdew
c and esat are mostly irrelevant, and don't cause any inconsistencies in
c logic. The biggest implication of this inconsistency is that we may
c not quite be conserving tdew when we take the diurnal cycle, even when
c we intend to do so. But I think the errors introduced by this
c inconsistency are small.
c
c Nevertheless, it is important to be aware of this inconsistency if you
c add any more equations converting between various forms of humidity.
c Here is an example of something that would NOT be correct: In this
c example, we try to use 100% RH at the coldest time of night rather
c than 99% RH as is currently done:
c
c   tdew_min = tmin(i)   ! BAD!!!!
c   ... other code here, including calculing tdew_now from tdew_min ...
c   qa(i) = qsat (e_from_tdew (tdew_now), psurf(i))
c
c When tdew_now = tdew_min, the above code will give different results
c for qa(i) than if you said qa(i) = qsat (esat (tmin(i))), which is
c what is intended here. The code can be corrected by changing the first
c line to:
c
c   tdew_min = tdew_from_e (esat(tmin(i)))
c 
c Then the errors introduced by the inconsistency between e_from_tdew
c and esat will cancel out. In the above line, tdew_min will differ
c slightly from tmin, but the differences will be corrected when we
c later take e_from_tdew (tdew_now).
c
      use compar
c
      implicit none
c
      include 'comsat.h'
c
c Arguments
c
      real tdew    ! dew point temperature (K)

c Use different coefficients depending on whether tdew is above freezing
c or below freezing. (But note that the function is continuous around
c tdew = 273.16.)
      if (tdew .ge. 273.16) then
        e_from_tdew = tdew_a * exp(tdew_b_w * (tdew - 273.16) /
     >   ((tdew - 273.16) + tdew_c_w))
      else
        e_from_tdew = tdew_a * exp(tdew_b_i * (tdew - 273.16) /
     >   ((tdew - 273.16) + tdew_c_i))
      end if

      return
      end


