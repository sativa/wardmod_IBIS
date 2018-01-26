c
c #####   #    #   #   #   ####      #     ####   #        ####    ####    #   #
c #    #  #    #    # #   #          #    #    #  #       #    #  #    #    # #
c #    #  ######     #     ####      #    #    #  #       #    #  #          #
c #####   #    #     #         #     #    #    #  #       #    #  #  ###     #
c #       #    #     #    #    #     #    #    #  #       #    #  #    #     #
c #       #    #     #     ####      #     ####   ######   ####    ####      #
c
c ---------------------------------------------------------------------
      subroutine stomata
c ---------------------------------------------------------------------
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsno
      use comsoi
      use comveg
      use com1d
      use comsum
      use comcrop
      use comnitr
      use compft
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      integer i, j, idc,
     >        iday,
     >        imonth,
     >        iyear
c
      real rwork,  ! 3.47e-03 - 1. / tu(i)
     >     tau,    ! 
     >     tleaf,  ! leaf temp in celcius
     >     tempvm, !
     >     zweight !

      real esat12, ! saturation vapor pressure in upper canopy air 
     >     e12,    ! vapor pressure in upper canopy air
     >     rh12,   ! relative humidity in upper canopy air 
     >     esat34, ! saturation vapor pressure in lower canopy air
     >     e34,    ! vapor pressure in lower canopy air 
     >     rh34,   ! relative humidity in lower canopy air 
     >     gbco2u, ! bound. lay. conductance for CO2 in upper canopy
     >     gbco2l, ! bound. lay. conductance for CO2 in lower canopy
     >     gscub,  ! 
     >     gscuc,  !
     >     gscls,  !
     >     gscl3,  !
     >     gscl4,  !
     >     gscc3,  !
     >     gscc4   !  
      real vmax, vmaxub, vmaxuc, vmaxls, vmaxl3, vmaxl4 
      real rdarkub, rdarkuc, rdarkls, rdarkl3, rdarkl4, rdarkc3, rdarkc4
      real agub, aguc, agls, agl3, agl4, agc3, agc4
      real anub, anuc, anls, anl3, anl4, anc3, anc4
      real duma, dumb, dumc, dume, dumq, dump
      real pxaiu, plaiu, pxail, plail
      real cscub, cscuc, cscls, cscl3, cscl4, cscc3, cscc4
      real extpar, scale
      real stressc3c, stressc4c
c
      real kc,     ! co2 kinetic parameter (mol/mol)
     >     ko,     ! o2  kinetic parameter (mol/mol)
*    >     ko15,   ! o2  kinetic parameter (mol/mol) at 15 degrees C
     >     kco2,   ! initial c4 co2 efficiency (mol-co2/m**2/s)
     >     je,     ! 'light limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jc,     ! 'rubisco limited' rate of photosynthesis (mol-co2/m**2/s)
     >     js,     ! 'sucrose limited' rate of photosynthesis (mol-co2/m**2/s)
     >     ji,     ! 'co2 limited' rate of photosynthesis (mol-co2/m**2/s)
     >     jp,     ! model-intermediate rate of photosynthesis (mol-co2/m**2/s)
     >     gamstar,! gamma*, the co2 compensation points for c3 plants
     >     q10,
     >     tkco2,
     >     absorb
c
c model parameters
c
c intrinsic quantum efficiency for c3 and c4 plants (dimensionless)
c
*      real alpha3, alpha4
c
*      data alpha3  /0.080/
*      data alpha4  /0.050/
c
c co2/o2 specificity ratio at 15 degrees C (dimensionless)
c
*      real tau15
c
*      data tau15 /4500.0/     
c
c o2/co2 kinetic parameters (mol/mol)
c
*      real kc15
c
c     data kc15 /1.5e-04/ 
c     data ko15 /2.5e-01/ 
c
c leaf respiration coefficients
c
*      real gammaub, gammauc, gammals, gammal3, gammal4,
*     >     gammac3, gammac4
c
*      data gammaub /0.0150/   ! broadleaf trees
*      data gammauc /0.0150/   ! conifer trees
*      data gammals /0.0150/   ! shrubs
*      data gammal3 /0.0150/   ! c3 grasses
*      data gammal4 /0.0300/   ! c4 grasses
*      data gammac3 /0.0150/   ! c3 crops
*      data gammac4 /0.0100/   ! c4 crops - corn  Amthor, 1984
c
c 'm' coefficients for stomatal conductance relationship
c
*      real coefmub, coefmuc, coefmls, coefml3, coefml4,
*     >     coefmc3, coefmc4
c
*     data coefmub /10.0/     ! broadleaf trees
*     data coefmuc / 6.0/     ! conifer trees
*     data coefmls / 9.0/     ! shrubs
*     data coefml3 / 9.0/     ! c3 grasses
*     data coefml4 / 4.0/     ! c4 grasses
*     data coefmc3 / 9.0/     ! c3 crops
*     data coefmc4 / 4.0/     ! c4 crops - corn (Collatz et al. 1992)
c
c 'b' coefficients for stomatal conductance relationship 
c (minimum conductance when net photosynthesis is zero)
c
*      real coefbub, coefbuc, coefbls, coefbl3, coefbl4,
*     >     coefbc3, coefbc4
c
*      data coefbub /0.010/    ! broadleaf trees
*      data coefbuc /0.010/    ! conifer trees
*      data coefbls /0.010/    ! shrubs
*      data coefbl3 /0.010/    ! c3 grasses
*      data coefbl4 /0.040/    ! c4 grasses
*      data coefbc3 /0.010/    ! c3 crops - soybean
c      data coefbc4 /0.080/    ! c4 crops - corn (Collatz et al. 1992)
*      data coefbc4 /0.030/    ! c4 crops - corn - Cupid model
c
c absolute minimum stomatal conductances
c
*      real gsubmin, gsucmin, gslsmin, gsl3min, gsl4min,
*     >     gsc3min, gsc4min
c
*      data gsubmin /0.00001/  ! broadleaf trees
*      data gsucmin /0.00001/  ! conifer trees
*      data gslsmin /0.00001/  ! shrubs
*      data gsl3min /0.00001/  ! c3 grasses
*      data gsl4min /0.00001/  ! c4 grasses
*      data gsc3min /0.00001/  ! c3 crops
*      data gsc4min /0.00001/  ! c4 crops
c
c photosynthesis coupling coefficients (dimensionless)
c
*      real theta3, theta4, beta4, thetac4, betac4,
*     >     thetac3, betac3
c
*      data theta3 /0.950/     ! c3 photosynthesis
*      data beta3  /0.990/     ! c3 photosynthesis
*      data theta4 /0.970/     ! c4 photosynthesis (not crops)
*      data beta4  /0.800/     ! c4 photosynthesis (not crops)
c
c photosynthesis coupling coefficients for c4 crops 
c from Collatz et al. 1992
c
c
*      data thetac4 /0.970/    ! c4 crop photosynthesis - maize 
*      data betac4  /0.800/    ! c4 crop photosynthesis - maize 
c
c photosynthesis coupling coefficients for c3 crops
c from Cupid model (Norman, 1983)
c
*      data thetac3 /0.950/    ! c3 crop photosynthesis - soybean 
*      data betac3  /0.990/    ! c3 crop photosynthesis - soybean
c
c
c maximum values for ci (for model stability)
c
*     real cimax
c
*     data cimax /2000.e-06/  ! maximum values for ci
c
c crop parameters
c
*      real
*     >     lotemp(npft),      ! low temperature threshold in tempvm equation
*     >     hitemp(npft),      ! high temperature threshold in tempvm equation 
*     >     drought(npft),     ! crop sensitivity to drought parameter 
*     >     f1(npft),          ! constant used in tempvm equations 
*     >     f2(npft)           ! constant used in tempvm equations
c   
c include water vapor functions
c
      include 'comsat.h'
c
c ---------------------------------------------------------------------
c * * * upper canopy physiology calculations * * *
c ---------------------------------------------------------------------
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c only perform calculations if crops are not planted
c
c       if (icropsum(i) .eq. 0.) then
c calculate physiological parameter values which are a function of temperature
c
        rwork = 3.47e-03 - 1. / tu(i)
c
        tau = tau15 * exp(-4500.0 * rwork)
        kc  = kc15  * exp( 6000.0 * rwork)
        ko  = ko15  * exp( 1500.0 * rwork)
c
        tleaf = tu(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * (  5.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c upper canopy gamma-star values (mol/mol)
c
        gamstar = o2conc / (2. * tau)
c
c constrain ci values to acceptable bounds -- to help ensure numerical stability
c
        ciub(i) = max (1.05 * gamstar, min (cimax, ciub(i)))
        ciuc(i) = max (1.05 * gamstar, min (cimax, ciuc(i)))
c
c calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
c
        gbco2u = min (10.0, max (0.1, su(i) * 25.5))
c 
c calculate the relative humidity in the canopy air space
c with a minimum value of 0.30 to avoid errors in the 
c physiological calculations
c
        esat12 = esat (t12(i))
        e12 = e_from_q (q12(i), psurf(i))
        rh12   = max (0.30, e12 / esat12)
c
c ---------------------------------------------------------------------
c broadleaf (evergreen & deciduous) tree physiology 
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c tropical broadleaf trees          60.0 e-06 mol/m**2/sec
c warm-temperate broadleaf trees    40.0 e-06 mol/m**2/sec
c temperate broadleaf trees         25.0 e-06 mol/m**2/sec
c boreal broadleaf trees            25.0 e-06 mol/m**2/sec
c
        if (exist(i,1).gt.0.5) then
          vmaxub = vmax_pft(1)     
        else if (exist(i,3).gt.0.5) then
          vmaxub = vmax_pft(3)     
        else 
          vmaxub = vmax_pft(5) 
        endif
c
c vmax and dark respiration for current conditions
c
        vmax  = vmaxub * tempvm * stresstu(i)
        rdarkub = gammaub * vmaxub * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparu(i) * 4.59e-06 * alpha3 * (ciub(i) - gamstar) / 
     >       (ciub(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (ciub(i) - gamstar) / 
     >       (ciub(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c       calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c
        js = vmax / 2.2
c
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agub = min (dumq/duma, dumc/dumq)
        anub = agub - rdarkub
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csub(i) = 0.5 * (csub(i) + co2conc - anub / gbco2u)
        csub(i) = max (1.05 * gamstar, csub(i))
c
c calculate new value of gs using implicit scheme
c
        gsub(i) = 0.5 * (gsub(i)  +  (coefmub * anub * rh12 / csub(i) + 
     >                                coefbub * stresstu(i)))
c
        gsub(i) = max (gsubmin, coefbub * stresstu(i), gsub(i))
c
c calculate new value of ci using implicit scheme
c
        ciub(i) = 0.5 * (ciub(i) + csub(i) - 1.6 * anub / gsub(i))
        ciub(i) = max (1.05 * gamstar, min (cimax, ciub(i)))
c
c ---------------------------------------------------------------------
c conifer tree physiology 
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
c temperate conifer trees           30.0 e-06 mol/m**2/sec
c boreal conifer trees              20.0 e-06 mol/m**2/sec
c
        if (exist(i,4).gt.0.5) then
          vmaxuc = vmax_pft(4)
        else 
          vmaxuc = vmax_pft(6)
        endif
c
c vmax and dark respiration for current conditions
c
        vmax  = vmaxuc * tempvm * stresstu(i)
        rdarkuc = gammauc * vmaxuc * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparu(i) * 4.59e-06 * alpha3 * (ciuc(i) - gamstar) / 
     >       (ciuc(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (ciuc(i) - gamstar) / 
     >       (ciuc(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c       
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c       
        js = vmax / 2.2
c 
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c       
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        aguc = min (dumq/duma, dumc/dumq) 
        anuc = aguc - rdarkuc
c
c       if (i .eq. 1 .and. anuc .gt. 0.0) write(*,*) anuc*1e+06 
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csuc(i) = 0.5 * (csuc(i) + co2conc - anuc / gbco2u)
        csuc(i) = max (1.05 * gamstar, csuc(i))
c
c calculate new value of gs using implicit scheme
c
        gsuc(i) = 0.5 * (gsuc(i)  +  (coefmuc * anuc * rh12 / csuc(i) + 
     >                                coefbuc * stresstu(i)))
c
        gsuc(i) = max (gsucmin, coefbuc * stresstu(i), gsuc(i))
c
c calculate new value of ci using implicit scheme
c
        ciuc(i) = 0.5 * (ciuc(i) + csuc(i) - 1.6 * anuc / gsuc(i))
        ciuc(i) = max (1.05 * gamstar, min (cimax, ciuc(i)))
c
c ---------------------------------------------------------------------
c upper canopy scaling
c ---------------------------------------------------------------------
c
c the canopy scaling algorithm assumes that the net photosynthesis
c is proportional to absored par (apar) during the daytime. during night,
c the respiration is scaled using a 10-day running-average daytime canopy
c scaling parameter.
c
c apar(x) = A exp(-k x) + B exp(-h x) + C exp(h x)
c an(x) is proportional to apar(x)
c
c therefore, an(x) = an(0) * apar(x) / apar(0)
c an(x) = an(0) * (A exp(-k x) + B exp(-h x) + C exp(h x)) / 
c                 (A + B + C)
c
c this equation is further simplified to
c an(x) = an(0) * exp (-extpar * x)
c
c an(0) is calculated for a sunlit leaf at the top of the canopy using
c the full-blown plant physiology model (Farquhar/Ball&Berry, Collatz).
c then the approximate par extinction coefficient (extpar) is calculated
c using parameters obtained from the two-stream radiation calculation.
c
c an,canopy avg.= integral (an(x), from 0 to xai) / lai
c               = an(0) * (1 - exp (-extpar * xai )) / (extpar * lai)
c
c the term '(1 - exp (-extpar * xai )) / lai)' scales photosynthesis from leaf
c to canopy level (canopy average) at day time. A 10-day running mean of this
c scaling parameter (weighted by light) is then used to scale the respiration
c during night time.
c
c once canopy average photosynthesis is calculated, then the canopy average
c stomatal conductance is calculated using the 'big leaf approach',i.e. 
c assuming that the canopy is a big leaf and applying the leaf-level stomatal
c conductance equations to the whole canopy.
c
c calculate the approximate par extinction coefficient:
c
c extpar = (k * A + h * B - h * C) / (A + B + C)
c
        extpar = (termu(i,6) * scalcoefu(i,1) +
     >            termu(i,7) * scalcoefu(i,2) -
     >            termu(i,7) * scalcoefu(i,3)) /
     >            max (scalcoefu(i,4), epsilon)
c
        extpar = max (1.e-1, min (1.e+1, extpar))
c
c calculate canopy average photosynthesis (per unit leaf area):
c
        pxaiu = extpar * (lai(i,2) + sai(i,2))
        plaiu = extpar *  lai(i,2)
c
c scale is the parameter that scales from leaf-level photosynthesis to
c canopy average photosynthesis
c CD : replaced 24 (hours) by 86400/dtime for use with other timestep
c
         zweight = exp(-1. / (10.0 * 86400. / dtime))
c
c for non-zero lai
c
        if (plaiu.gt.0.0) then
c
c day-time conditions, use current scaling coefficient
c
          if (topparu(i).gt.10.) then
c
            scale = (1. - exp(-pxaiu)) / plaiu
c
c update 10-day running mean of scale, weighted by light levels
c
            a10scalparamu(i) = zweight * a10scalparamu(i) + 
     >                         (1. - zweight) * scale * topparu(i)
c
            a10daylightu(i)  = zweight * a10daylightu(i) + 
     >                         (1. - zweight) * topparu(i)
c
c night-time conditions, use long-term day-time average scaling coefficient
c
          else
c
            scale = a10scalparamu(i) / a10daylightu(i)
c
          endif
c
c if no lai present
c
        else
c
          scale = 0.0
c
        endif
c
c perform scaling on all carbon fluxes from upper canopy
c
        agcub(i) = agub * scale
        agcuc(i) = aguc * scale
c
        ancub(i) = anub * scale
        ancuc(i) = anuc * scale
c
c calculate diagnostic canopy average surface co2 concentration 
c (big leaf approach)
c
        cscub = max (1.05 * gamstar, co2conc - ancub(i) / gbco2u)
        cscuc = max (1.05 * gamstar, co2conc - ancuc(i) / gbco2u)
c
c calculate diagnostic canopy average stomatal conductance (big leaf approach)
c
        gscub = coefmub * ancub(i) * rh12 / cscub +
     >          coefbub * stresstu(i)
c
        gscuc = coefmuc * ancuc(i) * rh12 / cscuc +
     >          coefbuc * stresstu(i)
c
        gscub = max (gsubmin, coefbub * stresstu(i), gscub)
        gscuc = max (gsucmin, coefbuc * stresstu(i), gscuc)
c
c calculate total canopy and boundary-layer total conductance for 
c water vapor diffusion
c
        rwork = 1. / su(i)
        dump  = 1. / 0.029
c
        totcondub(i) = 1. / (rwork + dump / gscub)
        totconduc(i) = 1. / (rwork + dump / gscuc)
c
c multiply canopy photosynthesis by wet fraction - this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1 - fwetu(i)
c
        agcub(i) = rwork * agcub(i)
        agcuc(i) = rwork * agcuc(i)
c
        ancub(i) = rwork * ancub(i)
        ancuc(i) = rwork * ancuc(i)
c
c      endif  ! crop existence
c
 100  continue
c
c ---------------------------------------------------------------------
c * * * lower canopy physiology calculations * * *
c ---------------------------------------------------------------------
c
      do 200 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c calculate physiological parameter values which are a function of temperature
c
        rwork = 3.47e-03 - 1. / tl(i)
c
        tau = tau15 * exp(-5000.0 * rwork)
        kc  = kc15  * exp( 6000.0 * rwork)
        ko  = ko15  * exp( 1400.0 * rwork)
c
        tleaf = tl(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * (  5.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c lower canopy gamma-star values (mol/mol)
c
        gamstar = o2conc / (2. * tau)
c
c constrain ci values to acceptable bounds -- to help ensure numerical stability
c
        cils(i) = max (1.05 * gamstar, min (cimax, cils(i)))
        cil3(i) = max (1.05 * gamstar, min (cimax, cil3(i)))
        cil4(i) = max (0.0           , min (cimax, cil4(i)))
c
c constrain ci value to acceptable bounds for crops
c
        cic3(i) = max (1.05 * gamstar, min (cimax, cic3(i)))
        cic4(i) = max (0.0           , min (cimax, cic4(i)))
c
c calculate boundary layer parameters (mol/m**2/s) = su / 0.029 * 1.35
c
        gbco2l = min (10.0, max (0.1, sl(i) * 25.5))
c 
c calculate the relative humidity in the canopy air space
c with a minimum value of 0.30 to avoid errors in the 
c physiological calculations
c
        esat34 = esat (t34(i))
        e34 = e_from_q (q34(i), psurf(i))
        rh34   = max (0.30, e34 / esat34)
c
c only perform calculations below if crops are not planted
c
c      if (icropsum(i) .eq. 0.) then
c
c ---------------------------------------------------------------------
c shrub physiology
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
        vmaxls = vmax_pft(9) 
c 
c vmax and dark respiration for current conditions
c
        vmax  = vmaxls * tempvm * stresstl(i)
        rdarkls = gammals * vmaxls * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha3 * (cils(i) - gamstar) / 
     >       (cils(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cils(i) - gamstar) / 
     >       (cils(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c       
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c       
        js = vmax / 2.2
c 
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c       
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agls = min (dumq/duma, dumc/dumq)
        anls = agls - rdarkls
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csls(i) = 0.5 * (csls(i) + co2conc - anls / gbco2l)
        csls(i) = max (1.05 * gamstar, csls(i))
c
c calculate new value of gs using implicit scheme
c
        gsls(i) = 0.5 * (gsls(i) + coefmls * anls * rh34 / csls(i) +
     >                             coefbls * stresstl(i))
c
        gsls(i) = max (gslsmin, coefbls * stresstl(i), gsls(i))
c
c calculate new value of ci using implicit scheme
c
        cils(i) = 0.5 * (cils(i) + csls(i) - 1.6 * anls / gsls(i))
        cils(i) = max (1.05 * gamstar, min (cimax, cils(i)))
c
c ---------------------------------------------------------------------
c c3 grass physiology
c ---------------------------------------------------------------------
c 
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
        vmaxl3 = vmax_pft(12)
c 
c vmax and dark respiration for current conditions
c
        vmax  = vmaxl3 * tempvm * stresstl(i)
        rdarkl3 = gammal3 * vmaxl3 * tempvm
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha3 * (cil3(i) - gamstar) / 
     >       (cil3(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cil3(i) - gamstar) / 
     >       (cil3(i) + kc * (1. + o2conc / ko))
c
c solution to quadratic equation
c
        duma = theta3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c       
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c       
        js = vmax / 2.2
c 
c solution to quadratic equation
c
        duma = beta3
        dumb = jp + js
        dumc = jp * js
c       
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agl3 = min (dumq/duma, dumc/dumq)
        anl3 = agl3 - rdarkl3
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csl3(i) = 0.5 * (csl3(i) + co2conc - anl3 / gbco2l)
        csl3(i) = max (1.05 * gamstar, csl3(i))
c
c calculate new value of gs using implicit scheme
c
        gsl3(i) = 0.5 * (gsl3(i) + coefml3 * anl3 * rh34 / csl3(i) +
     >                   coefbl3 * stresstl(i))
c
        gsl3(i) = max (gsl3min, coefbl3 * stresstl(i), gsl3(i))
c
c calculate new value of ci using implicit scheme
c
        cil3(i) = 0.5 * (cil3(i) + csl3(i) - 1.6 * anl3 / gsl3(i))
        cil3(i) = max (1.05 * gamstar, min (cimax, cil3(i)))
c
c ---------------------------------------------------------------------
c c4 grass physiology
c ---------------------------------------------------------------------
c
c nominal values for vmax of top leaf at 15 C (mol-co2/m**2/s)
c
        vmaxl4 = vmax_pft(11)
c
c calculate the parameter values which are a function of temperature
c
        rwork = 3.47e-03 - 1. / tl(i)
c
        tleaf = tl(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /
     >           ((1.0 + exp(0.40 * ( 10.0 - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - 50.0))))
c
c vmax and dark respiration for current conditions
c
        vmax  = vmaxl4 * tempvm * stresstl(i)
        rdarkl4 = gammal4 * vmaxl4 * tempvm
c
c initial c4 co2 efficiency (mol/m**2/s)
c
        kco2 = 18.0e+03 * vmax
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha4
c
c 'rubisco limited' rate of photosynthesis
c
        jc = vmax
c
c solve for intermediate photosynthesis rate
c
        duma = theta4
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
c
        ji = kco2 * cil4(i)
c
c solution to quadratic equation
c
        duma = beta4
        dumb = jp + ji
        dumc = jp * ji
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agl4 = min (dumq/duma, dumc/dumq)
        anl4 = agl4 - rdarkl4
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c CD: For numerical stability (to avoid division by zero in gsl4), 
c csl4 is limited to 1e-8 mol_co2/mol_air.
c  
        csl4(i) = 0.5 * (csl4(i) + co2conc - anl4 / gbco2l)
        csl4(i) = max (1.e-8, csl4(i))
c
c calculate new value of gs using implicit scheme
c
        gsl4(i) = 0.5 * (gsl4(i) + coefml4 * anl4 * rh34 / csl4(i) +
     >                   coefbl4 * stresstl(i))
c
        gsl4(i) = max (gsl4min, coefbl4 * stresstl(i), gsl4(i))
c
c calculate new value of ci using implicit scheme
c
        cil4(i) = 0.5 * (cil4(i) + csl4(i) - 1.6 * anl4 / gsl4(i))
        cil4(i) = max (0.0, min (cimax, cil4(i)))
c
c      endif ! crop existence
c
c --------------------------------------------------------------------
c identify crops that are planted 
c ---------------------------------------------------------------------
c
c temperature response curve constants which vary
c
*       f1(13) = 0.40
*       f1(14) = 0.40
*       f1(15) = 0.40
*       f2(13) = 0.40
*       f2(14) = 0.40
*       f2(15) = 0.40
c
c hi and low temperature thresholds (C)
c  
*       lotemp(13) = 5.0
*       lotemp(14) = 6.0
*       lotemp(15) = 0.0 
c
*       hitemp(13) = 40.0
*       hitemp(14) = 50.0
*       hitemp(15) = 38.0
c
c vmax values 15 C base
c
*       vmax_pft(13) = 65.0e-06  ! soybean vmax  Wullschleger, 1993
*       vmax_pft(14) = 70.0e-06  ! maize vmax
*       vmax_pft(15) = 60.0e-06  ! wheat vmax    Wullschleger, 1993 
c
c drought sensitivity - to account for differences between
c soybeans and other cereal crops (e.g., maize, wheat)
c
*       drought(13)  = 1.25
*       drought(14)  = 1.00
*       drought(15)  = 1.00
c       
       idc = 0
       do 250 j = scpft, ecpft 
c
c get index for current cropping practice - could have more than
c one crop existing in grid cell during the year for multiple
c cropping - but only one would be in live vegetation stage 
c
         if (exist(i,j) .eq. 1. .and. croplive(i,j) .gt. 0) then
            idc = j 
         endif
 250   continue 
c
c --------------------------------------------------------------------
c c3 crops physiology (soybean, wheat)
c ---------------------------------------------------------------------
       if (idc .eq. 13 .or. idc .eq. 15) then  ! soybean or wheat  
c
         rwork = 3.47e-03 - 1. / tl(i)
c
         tleaf = tl(i) - 273.16
c
         q10 = 2.0
c
c vmax and dark respiration for current conditions
c
         tempvm = exp(3500.0 * rwork ) /
     >            ((1.0 + exp(0.40 * (  lotemp(idc) - tleaf))) * 
     >            (1.0 + exp(0.40 * (tleaf - hitemp(idc)))))
c
c CJK 5/18/07        tempvm = q10**((tleaf-15.0)/10.0) /    ! Collatz approach
c     >           ((1.0 + exp(f1(idc) * (lotemp(idc) - tleaf))) * 
c     >            (1.0 + exp(f2(idc) * (tleaf - hitemp(idc)))))
c
c adjust drystress factor imposed on soybeans - on a scale of
c 0 (less) - 1.0 (more), these have a 0.8 rating compared to 0.65 for maize 
c and for wheat
c from Penning de Vries, "Simulation of ecophysiological processes of
c growth in several annual crops"
c
c make average stress factor 25% higher to account for difference 
c
c       stressc3c = min(1.0, stresstl(i) * drought(idc))
c       NEW CJK 9-24-04
c
c	stressc3c = min(1.0, drought(idc))
        stressc3c = min(1.0, stresstl(i))
        if (nstress .eq. 0) then  ! no N stress
          vmax      = max(0., vmax_pft(idc) * tempvm * 
     >                min(stressc3c, croplive(i,idc)))
        else  ! include N stress
          vmax      = max(0., vmax_pft(idc) * tempvm * 
     >                min(stressc3c, stressn(i,idc), croplive(i,idc)))
        end if
c       vmax      = max(0., vmax_pft(idc) * tempvm * croplive(i,idc))
c       vmax      = max(0., vmax_pft(idc) * tempvm * stressc3c * croplive(i,idc))
        rdarkc3   = gammac3 * vmax_pft(idc) * tempvm * croplive(i,idc)
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha3 * (cic3(i) - gamstar) / 
     >       (cic3(i) + 2. * gamstar)
c
c 'rubisco limited' rate of photosynthesis (mol/m**2/s)
c
        jc = vmax * (cic3(i) - gamstar) / 
     >       (cic3(i) + kc * (1. + o2conc / ko))
c
        duma = thetac3
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the intermediate photosynthesis rate (mol/m**2/s)
c
        jp = min (dumq/duma, dumc/dumq)
c       
c 'sucrose synthesis limited' rate of photosynthesis (mol/m**2/s)
c       
        js = vmax / 2.2
c 
c solution to quadratic equation
c
        duma = betac3
        dumb = jp + js
        dumc = jp * js
c       
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agc3 = min (dumq/duma, dumc/dumq)
        anc3 = agc3 - rdarkc3
c
c apply stress functions to net photosynthesis rate
c
c        anc3 = anc3 * stresstl(i)  ! CJK 6/20/2004  
c        open (23, file='physiology.dat', status='unknown')
c        if (i.eq.4) write(23,30) je*1e+06, vmax*1e+06, jc*1e+06, anc3*1e+06, tleaf
c30      format(5f8.2) 
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csc3(i) = 0.5 * (csc3(i) + co2conc - anc3 / gbco2l)
        csc3(i) = max (1.05 * gamstar, csc3(i))
c
c calculate new value of gs using implicit scheme
c
        gsc3(i) = 0.5 * (gsc3(i) + coefmc3 * anc3 * rh34 / csc3(i) +
     >                   coefbc3 * stressc3c)
c
        gsc3(i) = max (gsc3min, coefbc3 * stressc3c, gsc3(i))
c
c calculate new value of ci using implicit scheme
c
        cic3(i) = 0.5 * (cic3(i) + csc3(i) - 1.6 * anc3 / gsc3(i))
        cic3(i) = max (1.05 * gamstar, min (cimax, cic3(i)))
c
      else
        agc3    = 0.
        anc3    = 0.     
        csc3(i) = 0.
        gsc3(i) = 0.
        cic3(i) = 0.
c
      endif  ! c3 crops
c
c ---------------------------------------------------------------------
c c4 crop physiology (corn)
c
c modified by CJK to follow Collatz et al. (1992) parameterizations
c that were based on corn data in the field
c ---------------------------------------------------------------------
c
      if (idc .eq. 14) then   ! maize    
c
c calculate the parameter values which are a function of temperature
c
        q10   = 2.0
        rwork = 3.47e-03 - 1. / tl(i)
        tleaf = tl(i) - 273.16
c
        tempvm = exp(3500.0 * rwork ) /        ! original Agro-IBIS 1/26/2006
c  CJK 5/18/07       tempvm = q10**((tleaf-15.0)/10.0) /    ! Collatz approach
     >           ((1.0 + exp(f1(idc) * (lotemp(idc) - tleaf))) * 
     >            (1.0 + exp(f2(idc) * (tleaf - hitemp(idc)))))
c
c temperature effect on dark respiration - changed to 15 C base
c 
c         tresp  = q10**((tleaf - 15.0)/10.0) /
c     >            (1.0 + exp(1.30 * (tleaf - 55.0)))  
c
c vmax and dark respiration for current conditions
c add nitrogen stress to c4 crops 
c
c       stressc4c = min(1.0, stresstl(i) * drought(idc))
c       stressc4c = min(1.0, drought(idc))
        stressc4c = min(1.0, stresstl(i))   ! CJK 8/23/2006 
        if (nstress .eq. 0) then  ! no N stress
          vmax    = max(0., vmax_pft(idc) * tempvm * 
     >              min(stressc4c, croplive(i,idc)))
        else                    ! include N stress
          vmax    = max(0., vmax_pft(idc) * tempvm * 
     >              min(stressc4c, stressn(i,idc), croplive(i,idc)))
        end if
c       vmax    = max(0., vmax_pft(idc)   * tempvm * stresstl(i) * croplive(i,idc))
c
        rdarkc4 = gammac4 * vmax_pft(idc) * tempvm * croplive(i,idc)
c       rdarkc4 = gammac4 * vmax_pft(idc) * tempvm 
c
c from Collatz et al. (1992) - changed to 15 C base 
c equation 5B
c
        tkco2 = q10**((tleaf - 15.0)/10)
c
c initial c4 co2 efficiency (mol/m**2/s)
c
        kco2 = 18.0e+03 * tkco2 * vmax
c       kco2 = 18.0e+03 * vmax
c
c 'light limited' rate of photosynthesis (mol/m**2/s)
c
        je = topparl(i) * 4.59e-06 * alpha4 ! original C4 all plants - collatz 
c CJK 5/18/07       je = topparl(i) * 4.59e-06  * 0.067 ! needed to increase efficiency of corn plants 
c
c 'rubisco limited' rate of photosynthesis
c
        jc = vmax
c
c solve for intermediate photosynthesis rate
c
        duma = thetac4
        dumb = je + jc
        dumc = je * jc
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
        jp = min (dumq/duma, dumc/dumq)
c
c 'carbon dioxide limited' rate of photosynthesis (mol/m**2/s)
c
        ji = kco2 * cic4(i)
c
        duma = betac4
        dumb = jp + ji
        dumc = jp * ji
c
        dume = max (dumb**2 - 4. * duma * dumc, 0.)
        dumq = 0.5 * (dumb + sqrt(dume)) + 1.e-15
c
c calculate the net photosynthesis rate (mol/m**2/s)
c
        agc4 = min (dumq/duma, dumc/dumq)
        anc4 = agc4 - rdarkc4
c
c apply stress functions to net photosynthesis rate
c
c        anc4 = anc4 * max(0.0, stresstl(i))  ! CJK 6/20/2004 
c
c calculate co2 concentrations and stomatal condutance values
c using simple iterative procedure
c
c weight results with the previous iteration's values -- this
c improves convergence by avoiding flip-flop between diffusion
c into and out of the stomatal cavities
c
c calculate new value of cs using implicit scheme
c
        csc4(i) = 0.5 * (csc4(i) + co2conc - anc4 / gbco2l)
        csc4(i) = max (0.0, csc4(i))
c
c calculate new value of gs using implicit scheme
c
        gsc4(i) = 0.5 * (gsc4(i) + coefmc4 * anc4 * rh34 / csc4(i) +
     >                   coefbc4 * stressc4c)
c
        gsc4(i) = max (gsc4min, coefbc4 * stressc4c, gsc4(i))
c
c calculate new value of ci using implicit scheme
c
        cic4(i) = 0.5 * (cic4(i) + csc4(i) - 1.6 * anc4 / gsc4(i))
        cic4(i) = max (0.0, min (cimax, cic4(i)))
c
      else
        agc4    = 0.
        anc4    = 0.     
        csc4(i) = 0.
        gsc4(i) = 0.
        cic4(i) = 0.
c
      endif    ! c4 crops
c ---------------------------------------------------------------------
c lower canopy scaling
c ---------------------------------------------------------------------
c
c calculate the approximate extinction coefficient
c
        extpar = (terml(i,6) * scalcoefl(i,1) + 
     >            terml(i,7) * scalcoefl(i,2) -
     >            terml(i,7) * scalcoefl(i,3)) /
     >            max (scalcoefl(i,4), epsilon)
c
        extpar = max (1.e-1, min (1.e+1, extpar))
c
c calculate canopy average photosynthesis (per unit leaf area):
c
        pxail = extpar * (lai(i,1) + sai(i,1))
        plail = extpar *  lai(i,1)
c
c scale is the parameter that scales from leaf-level photosynthesis to
c canopy average photosynthesis
c CD : replaced 24 (hours) by 86400/dtime for use with other timestep
c
        zweight = exp(-1. / (10.0 * 86400. / dtime))
c
c for non-zero lai
c
        if (plail.gt.0.0) then
c
c day-time conditions, use current scaling coefficient
c
          if (topparl(i).gt.10.) then
c
            scale = (1. - exp(-pxail)) / plail
c
c update 10-day running mean of scale, weighted by light levels
c
            a10scalparaml(i) = zweight * a10scalparaml(i) + 
     >                         (1. - zweight) * scale * topparl(i)
c
            a10daylightl(i)  = zweight * a10daylightl(i) + 
     >                         (1. - zweight) * topparl(i)
c
c night-time conditions, use long-term day-time average scaling coefficient
c
          else
c
            scale = a10scalparaml(i) / a10daylightl(i)
c
          endif
c
c if no lai present
c
        else
c
          scale = 0.0
c
        endif
c
c perform scaling on all carbon fluxes from lower canopy
c
        agcls(i) = agls * scale
        agcl4(i) = agl4 * scale
        agcl3(i) = agl3 * scale
        agcc3(i) = agc3 * scale
        agcc4(i) = agc4 * scale
c
        ancls(i) = anls * scale
        ancl4(i) = anl4 * scale
        ancl3(i) = anl3 * scale
        ancc3(i) = anc3 * scale
        ancc4(i) = anc4 * scale
c
c calculate canopy average surface co2 concentration
c CD: For numerical stability (to avoid division by zero in gscl4),
c cscl4 is limited to 1e-8 mol_co2/mol_air.
c
        cscls = max (1.05 * gamstar, co2conc - ancls(i) / gbco2l)
        cscl3 = max (1.05 * gamstar, co2conc - ancl3(i) / gbco2l)
        cscl4 = max (1.0e-08       , co2conc - ancl4(i) / gbco2l)
c
        if (idc .eq. 13 .or. idc .eq. 15) then
          cscc3 = max (1.05 * gamstar, co2conc - ancc3(i) / gbco2l)
          cscc4 = 0.
        else if (idc .eq. 14) then
          cscc4 = max (1.0e-08       , co2conc - ancc4(i) / gbco2l)
          cscc3 = 0.
        endif
c
c calculate canopy average stomatal conductance
c
        gscls = coefmls * ancls(i) * rh34 / cscls +
     >          coefbls * stresstl(i)
c
        gscl3 = coefml3 * ancl3(i) * rh34 / cscl3 +
     >          coefbl3 * stresstl(i)
c
        gscl4 = coefml4 * ancl4(i) * rh34 / cscl4 +
     >          coefbl4 * stresstl(i)
c
c c3/c4 crop physiology
c
        if (idc .eq. 13 .or. idc .eq. 15) then
          gscc3 = coefmc3 * ancc3(i) * rh34 / cscc3 +
     >            coefbc3 * stressc3c
          gscc4 = 0.
c
        else if (idc .eq. 14) then
          gscc4 = coefmc4 * ancc4(i) * rh34 / cscc4 +
     >            coefbc4 * stressc4c
          gscc3 = 0.  
        else
          gscc3 = 0.
          gscc4 = 0.
        endif
c
        gscls = max (gslsmin, coefbls * stresstl(i), gscls)
        gscl3 = max (gsl3min, coefbl3 * stresstl(i), gscl3)
        gscl4 = max (gsl4min, coefbl4 * stresstl(i), gscl4)
c
c c3/c4 crop physiology
c
        if (idc .eq. 13 .or. idc .eq. 15) then
          gscc3 = max (gsc3min, coefbc3 * stressc3c, gscc3)
          gscc4 = 0.
        else if (idc .eq. 14) then
          gscc4 = max (gsc4min, coefbc4 * stressc4c, gscc4)
          gscc3 = 0.
        else
          gscc3 = 0.
          gscc4 = 0.
        endif
c
c The following adjusts the above calculated values of ancl3, ancl4,
c agcl3, agcl4, gscl3, and gscl4 according to what percentage of the
c lower canopy is green by weighting the above calculations by greenfrac
c terms. Only the green portion of the canopy performs photosynthesis.
c Shrubs that have leaves have only green leaves since they are allowed
c to drop their leaves in the fall. C3 and C4 grasses may be either green
c or brown so they are the only terms that are adjusted.
c
c Scale value of ancl3, ancl4, gscl3, and gscl4 according to what fraction
c of the canopy is green
c
        ancl3(i) = ancl3(i) * greenfracl3(i)
        ancl4(i) = ancl4(i) * greenfracl4(i)
c
        agcl3(i) = agcl3(i) * greenfracl3(i)
        agcl4(i) = agcl4(i) * greenfracl4(i)
c
        gscl3 = gscl3 * greenfracl3(i)
        gscl4 = gscl4 * greenfracl4(i)
c
c Do similar adjustments for crops
c
        if (idc .eq. 13 .or. idc .eq. 15) then
          ancc3(i) = ancc3(i) * grnfraccrop(i,idc)
          agcc3(i) = agcc3(i) * grnfraccrop(i,idc)
          gscc3 = gscc3 * grnfraccrop(i,idc)
        else if (idc .eq. 14) then
          ancc4(i) = ancc4(i) * grnfraccrop(i,idc)
          agcc4(i) = agcc4(i) * grnfraccrop(i,idc)
          gscc4 = gscc4 * grnfraccrop(i,idc)
        end if          
c          
c
c calculate canopy and boundary-layer total conductance for water vapor diffusion
c
        rwork = 1. / sl(i)
        dump =  1. / 0.029
c
        totcondls(i) = 1. / (rwork + dump / gscls)
c
c Make sure that the calculation does not divide by zero if gscl3 or
c gscl4 are equal to zero
c
        if (gscl3 .gt. 0) then
          totcondl3(i) = 1. / ( rwork + dump / gscl3 )
        else
          totcondl3(i) = 0
        endif
c
        if (gscl4 .gt. 0) then
          totcondl4(i) = 1. / ( rwork + dump / gscl4 )
        else
          totcondl4(i) = 0
        endif
c
        if (gscc3 .gt. 0) then
          totcondc3(i) = 1. / ( rwork + dump / gscc3 )
        else
          totcondc3(i) = 0
        endif
c
        if (gscc4 .gt. 0) then
          totcondc4(i) = 1. / ( rwork + dump / gscc4 )
        else
          totcondc4(i) = 0
        endif
c
c
c multiply canopy photosynthesis by wet fraction -- this calculation is
c done here and not earlier to avoid using within canopy conductance
c
        rwork = 1. - fwetl(i)
c
        agcls(i) = rwork * agcls(i)
        agcl3(i) = rwork * agcl3(i)
        agcl4(i) = rwork * agcl4(i)
        agcc3(i) = rwork * agcc3(i)
        agcc4(i) = rwork * agcc4(i)
c
        ancls(i) = rwork * ancls(i)
        ancl3(i) = rwork * ancl3(i)
        ancl4(i) = rwork * ancl4(i)
        ancc3(i) = rwork * ancc3(i)
        ancc4(i) = rwork * ancc4(i)
c
 200  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine drystress
c ---------------------------------------------------------------------
c
c modified by CJK 1/25/2005 for dynamic root water uptake adjustment
c approach of K. Li
c
c uses:
c
      use comgrid
      use compar
      use comsoi
      use comveg
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      integer i, k,   ! loop indices
     >        klay    ! layer of soil that is one third of total depth 
c
      real stressfac, ! to calculate moisture stress factor 
     >     awc,       ! available water content (fraction)
     >     znorm,     ! normalizing factor
     >     zwilt,     ! function of awc, =1 if awc = 1 (no stress)
     >     f1,
     >     f2,
     >     foptawc,
     >     afact,
     >     lambda,    ! exponent in root water uptake equation
     >     sdep,      ! total soil depth calculation
     >     tthird,    ! top third of soil depth
     >     botawc,
     >     topawc,
     >     botmxw,
     >     topmxw,
     >     topw,
     >     botw,
     >     mxwrate
c
c stressfac determines the 'strength' of the soil moisture
c stress on physiological processes
c
c strictly speaking, stresst* is multiplied to net photosynthesis 
c parameters used in the photosynthesis calculations
c
c stressfac determines the shape of the soil moisture response
c
c      stressfac = -5.0
c
c      znorm = 1.0 - exp(stressfac)
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c initialize stress parameter
c
        stresstl(i) = 0.0
        stresstu(i) = 0.0
c
c set values of lambda
c
c        sdep = 0.    ! total soil depth
c        do 105 k = 1, nsoilay
c          sdep = sdep + hsoi(k) 
c 105    continue 
c
c        tthird = sdep * 0.3333   ! calculate 1/3 of total depth
c        sdep = 0.
c        klay = 0
c        do 110 k = 1, nsoilay
c          sdep = sdep + hsoi(k)
c          if (sdep .le. tthird) klay = k
c 110    continue      
c
c        topw   = 0.
c        botw   = 0.
c        topmxw = 0.
c        botmxw = 0.
c        topawc = 0.
c        botawc = 0.
c        do 120 k = 1, klay 
c          topw   = topw + (wsoi(i,k)*(1-wisoi(i,k)) - swilt(i,k)) 
c          topmxw = topmxw + (sfield(i,k) - swilt(i,k)) 
c 120    continue
c        topawc = topw / topmxw 
c
c        do 130 k = klay+1, nsoilay 
c          botw   = botw + (wsoi(i,k)*(1-wisoi(i,k)) - swilt(i,k)) 
c          botmxw = botmxw + (sfield(i,k) - swilt(i,k)) 
c 130    continue
c        botawc = botw / botmxw 
c      
c assign lambda values
c
c        if (topawc .lt. 0.20 .and. botawc .gt. 0.5) then
c          lambda = 0.50 
c        else if (topawc .lt. 0.20 .and.
c     >           (botawc .ge. 0.2 .and. botawc .le. 0.5)) then
c          lambda = 0.75
c        else if ((topawc .ge. 0.2 .and. topawc .le. 0.5) .and.
c     >            botawc .gt. 0.5) then
c          lambda = 0.75         
c        else if ((topawc .ge. 0.2 .and. topawc .le. 0.5) .and.
c     >            botawc .lt. 0.2) then
c          lambda = 1.25         
c        else if (topawc .gt. 0.5 .and.
c     >           ( botawc .ge. 0.2 .and. botawc .le. 0.5)) then
c          lambda = 1.25         
c        else if (topawc .gt. 0.5 .and. botawc .lt. 0.2) then
c          lambda = 1.50         
c        else
c          lambda = 1.00
c        endif
c
c fraction of soil water uptake in each layer
c
        do 200 k = 1, nsoilay
c
          awc = min (1.0, max (0.0,
     >              (wsoi(i,k)*(1 - wisoi(i,k))   - swilt(i,k)) /
     >              (sfield(i,k) - swilt(i,k))
     >              )         )
c
c original IBIS formulation for zwilt 
c
c         zwilt = (1. - exp(stressfac * awc)) / znorm
c
c 
c J. Norman soil water availability factor (zwilt) 
c 1/25/2004
c
          zwilt = 1.0 - (alog(1+799.0*exp(-12.0*awc))/alog(800.))  
c
c update for each layer
c
           stressl(i,k) = froot(k,1) * max (0.0, min (1.0, zwilt))
           stressu(i,k) = froot(k,2) * max (0.0, min (1.0, zwilt))
c
c          stressl(i,k) = froot(k,1)**lambda * max (0.0, min (1.0, zwilt))
c          stressu(i,k) = froot(k,2)**lambda * max (0.0, min (1.0, zwilt))
c
c calculate maximum dimensionless water uptake rate (ranging from 0-1)
c
c          mxwrate = 1. - (1 + 1.3 * awc)**(-bex(i,k)) 
c          stressl(i,k) = stressl(i,k) * mxwrate
c          stressu(i,k) = stressu(i,k) * mxwrate
c
c integral over rooting profile
c
          stresstl(i) = stresstl(i) + stressl(i,k)
          stresstu(i) = stresstu(i) + stressu(i,k)
c
 200    continue
          stresstl(i) = min(1.0, max(0.001, stresstl(i)))
c
 100  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine co2 (co2init, co2conc, iyear)
c ---------------------------------------------------------------------
c
c     use comveg
c
      implicit none
c
c Arguments 
c
      integer iyear,   ! current year
     >        iyr
c
      real co2init,    ! input atmospheric co2 concentration
     >     co2conc     ! output " for year iyear   
c
c calculate co2 concentration for this year
c
      if (iyear.lt.1860) then
c
        co2conc = co2init
c
      else
c
c 1992 IPCC estimates
c
c        iyr = iyear - 1860 + 1
c        co2conc = (297.12 - 0.26716 * iyr +
c     >                      0.0015368 * iyr**2 +
c     >                      3.451e-5 * iyr**3) * 1.e-6
c
c 1996 IPCC estimates
c
       iyr = iyear - 1860 + 1
       co2conc = (303.514 - 0.57881 * iyr +
     >                      0.00622 * iyr**2 +
     >                      1.3e-5 * iyr**3) * 1.e-6
c
      end if
c
      return
      end
c
