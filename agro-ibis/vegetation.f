c
c #    #  ######   ####   ######   #####    ##     #####     #     ####   #    #
c #    #  #       #    #  #          #     #  #      #       #    #    #  ##   #
c #    #  #####   #       #####      #    #    #     #       #    #    #  # #  #
c #    #  #       #  ###  #          #    ######     #       #    #    #  #  # #
c  #  #   #       #    #  #          #    #    #     #       #    #    #  #   ##
c   ##    ######   ####   ######     #    #    #     #       #     ####   #    #
c
c ---------------------------------------------------------------------
      subroutine pheno(jday)
c ---------------------------------------------------------------------
c
c uses:
c
      use comgrid
      use compar
      use comatm
      use comsoi
      use comsum
      use comveg
      use comcrop
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      integer
     >  i,j,          !
     >  jday          !
c
      real
     >  ddays,        !
     >  ddfac,        !
     >  avglaiu,      ! average lai of upper canopy 
     >  avglail,      ! average lai of lower canopy 
     >  sthreshold,   !  
     >  sumcrop
c
c define 'drop days' -- number of days to affect phenology change
c
      ddays = 15.0
      ddfac = 1.0 / ddays
c
c soil temperature threshold for budburst
c soil temperature threshold is assumed to be 0 degrees C
c
      sthreshold = 273.16
c
c initialize this year's values of stsumu, stsuml, precipsum
c
      if (jday.eq.1) then
c
        call const (stsumu, npoi, 0.)
        call const (stsuml, npoi, 0.)
        call const (precipsum, npoi, 0.)
        call logicf (onflagu, npoi)
        call logicf (onflagl3, npoi)
        call logicf (onflagl4, npoi)
        call logict (offflagu, npoi)
        call logict (offflagl3, npoi)
        call logict (offflagl4, npoi)
c
      endif
c
c begin global grid
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c only do if crops are not growing
c
        if (icropsum(i) .eq. 0.0) then
c
c ---------------------------------------------------------------------
c * * * upper canopy onset phenology * * *
c ---------------------------------------------------------------------
c
        if ((offflagu(i)) .AND. (jday .lt. 243)) then
c
c sumcom is the threshold for upper canopy onset
c
          sumcom(i) = EXP(4.395 + 0.129 * Tavgann(i))
c
c determine if soil temperature summation is initiated
c if so, calculate onset summation stsumu and stsuml
c
          if (a11soiltd(i).lt.sthreshold) then
            stsumu(i)  = stsumu(i)
          else
            stsumu(i) = stsumu(i) + a11soiltd(i) - sthreshold
          endif
c
c determine if onset has occured
c
          if (stsumu(i).ge.sumcom(i)) then
            onflagu(i) = .TRUE.
            offflagu(i) = .FALSE.
          endif
c
        endif
c
c if onset has occured then determine leaf display
c
        if (onflagu(i)) then
          tempu(i) = min (1., tempu(i) + ddfac)
        endif
c
c ---------------------------------------------------------------------
c * * * upper canopy offset phenology * * *
c ---------------------------------------------------------------------
c
        if (onflagu(i)) then
          if (jday .ge. 243) then
            if (((daylength(i).le.685.) .AND. (a11soiltd(i).le.(17.15+273.16)))
     >        .OR. (a11soiltd(i).le.(2.+273.16))) then
c
              offflagu(i) = .TRUE.
              onflagu(i) = .FALSE.
c
            endif
          endif
        endif
c
c if offset has occured then determine leaf display
c
        if (offflagu(i)) then
          tempu(i) = max (0., tempu(i) - ddfac)
        endif
c
c ---------------------------------------------------------------------
c * * * lower canopy onset phenology * * *
c ---------------------------------------------------------------------
c
c C3 and C4 grasses in the lower canopy act a little different than
c the trees found above. These grasses keep their leaves all year so
c their LAI is constant and they do not have a temp factor in the following
c algorithm. Instead, leaves are either green (living during the growing
c season) or brown (dead during the winter). Greenness is determined by
c the terms greenfracl3, for C3 grasses, or greenfracl4, for C4 grasses.
c The greenfrac terms act similarly to the temp factors. When onset occurs,
c greenfrac increases from 0 to 1 (all brown to all green leaves). When
c offset occurs, greenfrac decreases from 1 to 0.
c Greenfracl3 and greenfracl4 are then passed to stomata (physiology.f)
c and twoset (radiation.f) where they affect photosynthesis and leaf optical
c properties, respectively
c
        if ((offflagl3(i) .OR. offflagl4(i)) .AND. (jday .lt. 243)) then
c
c accumulate precipitation summation
c
          precipsum(i) = precipsum(i) + precip(i)
          PPTsumcrit(i) = 0.15 * PPTavgann(i)
c
c determine if soil temperature summation is initiated
c if so, calculate onset summation
c
          if (a11soiltd(i).lt.sthreshold) then
            stsuml(i)  = stsuml(i)
          else
            stsuml(i) = stsuml(i) + a11soiltd(i) - sthreshold
          endif
c
c
c for warm grasses, summation threshold is 1320
c for cool grasses, summation threshold is 420
c
          if (precipsum(i) .ge. PPTsumcrit(i)) then
            if (stsuml(i) .ge. 420.) then
              onflagl3(i) = .TRUE.
              offflagl3(i) = .FALSE.
            endif
            if (stsuml(i) .ge. 1320.) then
              onflagl4(i) = .TRUE.
              offflagl4(i) = .FALSE.
            endif
          endif
        endif
c
c if onset has occured then determine leaf color
c templs is retained so that deciduous shrubs may lose their leaves
c
        if (onflagl3(i)) then
           greenfracl3(i) = min (1., greenfracl3(i) + ddfac)
        endif
        if (onflagl4(i)) then
           greenfracl4(i) = min (1., greenfracl4(i) + ddfac)
        endif
c
c ---------------------------------------------------------------------
c * * * lower canopy offset phenology * * *
c ---------------------------------------------------------------------
c
c This is tuned White et al. version that looks at the stress of the plants
c and determines if 'cold' conditions are met
c
        if (onflagl3(i) .OR. onflagl4(i)) then
          if (jday .ge. 243) then
            if ((stresstl(i) .lt. 0.27) .OR.
     >        ((a3tdmin(i)-273.16).le.tminavgann(i))) then
c
              offflagl3(i) = .TRUE.
              offflagl4(i) = .TRUE.
              onflagl3(i) = .FALSE.
              onflagl4(i) = .FALSE.
c
            endif
          endif
        endif
c
c if offset has occured then determine leaf display
c
        if (offflagl3(i)) then
           greenfracl3(i) = max (0., greenfracl3(i) - ddfac)
        endif
c
        if (offflagl4(i)) then
           greenfracl4(i) = max (0., greenfracl4(i) - ddfac)
        endif
c
c ---------------------------------------------------------------------
c * * * update lai and canopy fractions * * *
c ---------------------------------------------------------------------
c
c Here the leaf display of shrubs, templs, is set equal to that of
c trees, tempu, since it was determined that even though shrubs are
c in the lower canopy, their leaf display follows more closely that
c of trees than grasses.
c
        templs(i) = tempu(i)
c
c upper canopy single sided leaf area index (area-weighted)
c
        avglaiu = plai(i,1)             +
     >            plai(i,2)             +
     >            plai(i,3)             +
     >            plai(i,4)             +
     >            plai(i,5) * tempu(i)  +
     >            plai(i,6)             +
     >            plai(i,7) * tempu(i)  +
     >            plai(i,8) * tempu(i)
c
c upper canopy fractions
c
        frac(i,1) = plai(i,1)            / max (avglaiu, epsilon)
        frac(i,2) = plai(i,2)            / max (avglaiu, epsilon)
        frac(i,3) = plai(i,3)            / max (avglaiu, epsilon)
        frac(i,4) = plai(i,4)            / max (avglaiu, epsilon)
        frac(i,5) = plai(i,5) * tempu(i) / max (avglaiu, epsilon)
        frac(i,6) = plai(i,6)            / max (avglaiu, epsilon)
        frac(i,7) = plai(i,7) * tempu(i) / max (avglaiu, epsilon)
        frac(i,8) = plai(i,8) * tempu(i) / max (avglaiu, epsilon)
c
c lower canopy single sided leaf area index (area-weighted)
c
        avglail = plai(i,9)                  +
     >            plai(i,10) *     templs(i) +
     >            plai(i,11)                 +
     >            plai(i,12)
c
c lower canopy fractions
c templs is included in frac(i,10) to allow
c deciduous shrubs to drop their leaves
c All other pfts keep leaves year-round
c
        frac(i,9)  = plai(i,9)                        /
     >               max (avglail, epsilon)
c
        frac(i,10) = plai(i,10) * templs(i)           /
     >               max (avglail, epsilon)
c
        frac(i,11) = plai(i,11)                        /
     >               max (avglail, epsilon)
c
        frac(i,12) = plai(i,12)                        /
     >               max (avglail, epsilon)
c
c Find an average fraction of green vegetation in the lower canopy
c to be used in stomata and twoset subroutines
c
        greenfracl(i) = frac(i,9) + frac(i,10)                  +
     >                              greenfracl4(i) * frac(i,11) +
     >                              greenfracl3(i) * frac(i,12)
c
c calculate the canopy leaf area index using the fractional vegetation cover
c
        lai(i,1) = avglail / fl(i)
        lai(i,2) = avglaiu / fu(i)
c
c put a fix on canopy lais to avoid problems in physics
c
        lai(i,1) = min (lai(i,1), 12.0)
        lai(i,2) = min (lai(i,2), 12.0)
c
c ---------------------------------------------------------------------
c * * * update canopy height parameters * * *
c ---------------------------------------------------------------------
c
c update lower canopy height parameters
c
c note that they are based on vegetation fraction and not
c averaged over the entire gridcell
c
c        zbot(i,1)   =  0.05
c        ztop(i,1)   =  max (0.25, lai(i,1) * 0.25)
c        
c constrain ztop to be at least 0.5 meter lower than 
c zbot for upper canopy
c
c        ztop(i,1) = min (ztop(i,1), zbot(i,2) - 0.5)
c
c
      endif ! crop existence check 
c end of loop
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
      subroutine dynaveg (iyear, isimfire)
c ---------------------------------------------------------------------
c
      use comgrid
      use compar
      use comsoi
      use comsum
      use comveg
      use compft
      use comcrop
      use combcs,only: cmask
c
      implicit none
c
c Arguments
c
      integer isimfire,   ! fire switch
     >        iyear       ! absolute year (1995,eg)

c
c local variables
c
      integer
     >  i, j              ! gridcell counter
c
      real
     >  fireresist(npft), ! pft fire resistances
     >  fracupper,        !fraction of upper canopy burned by fire
     >  fraclower,        !fraction of lower canopy burnder by fire
     >  sapspeed,         ! in mm/day
     >  trans,            ! (2.5 mm/day) 
     >  saparea,          ! in m**2
     >  sapvolume,        ! in m**3
     >  denswood,         ! kg/m**3
     >  wood,             ! total amount of woody biomass in gridcell
     >  taufin            !
*    >  xminlai           !
c
c      real
c     >  tauleaf(npft), ! turnover time of carbon in leaves (years)
c     >  tauroot(npft), ! turnover time of carbon in fine roots (years)
c     >  tauwood(npft), ! turnover time of carbon in wood (years)
c     >  tauwood0(npft) ! normal (unstressed) turnover time
c

c ibis uses a small number of plant functional types:
c
c  1: tropical broadleaf evergreen tree
c  2: tropical broadleaf drought-deciduous trees
c  3: warm-temperate broadleaf evergreen tree
c  4: temperate conifer evergreen tree
c  5: temperate broadleaf cold-deciduous tree
c  6: boreal conifer evergreen tree
c  7: boreal broadleaf cold-deciduous tree
c  8: boreal conifer cold-deciduous tree
c  9: evergreen shrub
c 10: deciduous shrub
c 11: warm (c4) grass
c 12: cool (c3) grass
c 13: soybeans
c 14: maize 
c 15: spring and winter wheat
c
c ---------------------------------------------------------------------
c * * * specify biomass turnover parameters (years) * * *
c ---------------------------------------------------------------------
c
*      data tauleaf / 1.01,   ! tropical broadleaf evergreen trees
*     >               1.00,   ! tropical broadleaf drought-deciduous trees
*     >               1.00,   ! warm-temperate broadleaf evergreen trees
*     >               2.00,   ! temperate conifer evergreen trees
*     >               1.00,   ! temperate broadleaf cold-deciduous trees
*     >               2.50,   ! boreal conifer evergreen trees
*     >               1.00,   ! boreal broadleaf cold-deciduous trees
*     >               1.00,   ! boreal conifer cold-deciduous trees
*     >               1.50,   ! evergreen shrubs
*     >               1.00,   ! deciduous shrubs
*     >               1.25,   ! warm (c4) grasses
*     >               1.50,   ! cool (c3) grasses
*     >               999.0,  ! soybean
*     >               999.0,  ! maize 
*     >               999.0 / ! wheat
c
*      data tauwood0 / 25.0,  ! tropical broadleaf evergreen trees
*     >                25.0,  ! tropical broadleaf drought-deciduous trees
*     >                25.0,  ! warm-temperate broadleaf evergreen trees
*     >                50.0,  ! temperate conifer evergreen trees
*     >                50.0,  ! temperate broadleaf cold-deciduous trees
*     >               100.0,  ! boreal conifer evergreen trees
*     >               100.0,  ! boreal broadleaf cold-deciduous trees
*     >               100.0,  ! boreal conifer cold-deciduous trees
*     >                 5.0,  ! evergreen shrubs
*     >                 5.0,  ! deciduous shrubs
*     >               999.0,  ! warm (c4) grasses
*     >               999.0,  ! cool (c3) grasses
*     >               999.0,  ! soybean
*     >               999.0,  ! maize 
*     >               999.0 / ! wheat 
c
c begin global grid
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
        if (icropsum(i) .eq. 0.0) then
c
c ---------------------------------------------------------------------
c * * * initialize vegetation dynamics pools * * *
c ---------------------------------------------------------------------
c
c zero out litter fall fields
c
        falll(i) = 0.0
        fallr(i) = 0.0
        fallw(i) = 0.0
c
c zero out carbon lost due to disturbance
c 
        cdisturb(i) = 0.0
c
        wood = 0.001
c
c ---------------------------------------------------------------------
c * * * update npp, and pool losses  * * *
c ---------------------------------------------------------------------
c
c go through all the pfts
c
        do 110 j = 1, npft
c
c apply this year's existence arrays to npp
c
          aynpp(i,j)  = exist(i,j) * aynpp(i,j)
c
c determine above-ground npp for each plant type
c
          ayanpp(i,j) = (aleaf(i,j) + awood(i,j)) * aynpp(i,j)
c
c determine turnover rates for woody biomass:
c
c if pft can exist,    then tauwood = tauwood0 (normal turnover),
c if pft cannot exist, then tauwood = taufin years (to kill off trees)
c
c          taufin     = 5.0
           taufin     = tauwood0(j)/2.0
c
          tauwood(j) = tauwood0(j) - (tauwood0(j) - taufin) *
     >                               (1.0 - exist(i,j))
c
c assume a constant fine root turnover time
c
          tauroot(j) = 1.0
c
c determine litter fall rates
c
          falll(i) = falll(i) + cbiol(i,j) / tauleaf(j)
          fallr(i) = fallr(i) + cbior(i,j) / tauroot(j)
          fallw(i) = fallw(i) + cbiow(i,j) / tauwood(j)
c
c ---------------------------------------------------------------------
c * * * update biomass pools  * * *
c ---------------------------------------------------------------------
c
c update carbon reservoirs using an analytical solution
c to the original carbon balance differential equation
c
          cbiol(i,j) = cbiol(i,j) * exp(-1./tauleaf(j))  +
     >                 aleaf(i,j) * tauleaf(j) * max (0., aynpp(i,j)) *
     >                 (1. - exp(-1./tauleaf(j)))
c
          cbiow(i,j) = cbiow(i,j) * exp(-1./tauwood(j))  +
     >                 awood(i,j) * tauwood(j) * max (0., aynpp(i,j)) *
     >                 (1. - exp(-1./tauwood(j)))
c
          cbior(i,j) = cbior(i,j) * exp(-1./tauroot(j))  +
     >                 aroot(i,j) * tauroot(j) * max (0., aynpp(i,j)) *
     >                 (1. - exp(-1./tauroot(j)))

c
          if (j.le.8) wood = wood + max (0.0, cbiow(i,j))
c
 110    continue
c
c ---------------------------------------------------------------------
c * * * apply disturbances * * *
c ---------------------------------------------------------------------
c
c set fixed disturbance regime
c
        firefrac(i) = 0.00
c
c if isimfire true, call fire subroutine
c

        if (isimfire.eq.1) call fire(iyear) 

c Define pft fire resistances MM 4/16/2010. Taken from Thonicke 2004.

        fireresist(1) = 0.12
        fireresist(2) = 0.5
        fireresist(3) = 0.5
        fireresist(4) = 0.12
        fireresist(5) = 0.12
        fireresist(6) = 0.12
        fireresist(7) = 0.12
        fireresist(8) = 0.12
        fireresist(9) = 1.0
        fireresist(10) = 1.0
        fireresist(11) = 1.0
        fireresist(12) = 1.0

        do 120 j = 1, npft

c
c calculate biomass (vegetations) carbon lost to atmosphere   
c used to balance net ecosystem exchange  
c
          cdisturb(i) = cdisturb(i) + 
     >             cbiol(i,j) * (firefrac(i) * fireresist(j)) +
     >             cbiow(i,j) * (firefrac(i) * fireresist(j)) +
     >             cbior(i,j) * (firefrac(i) * fireresist(j))

          
c adjust biomass pools due to disturbances
c
          cbiol(i,j) = cbiol(i,j) * (1. - firefrac(i)*fireresist(j))
          cbiow(i,j) = cbiow(i,j) * (1. - firefrac(i)*fireresist(j))
          cbior(i,j) = cbior(i,j) * (1. - firefrac(i)*fireresist(j))

c constrain biomass fields to be positive
c
          cbiol(i,j) = max (0.0, cbiol(i,j))
          cbiow(i,j) = max (0.0, cbiow(i,j))
          cbior(i,j) = max (0.0, cbior(i,j))
c
c maintain minimum value of leaf carbon in areas that plants exist
c
c         xminlai = 0.010
c
          cbiol(i,j) = max (exist(i,j) * xminlai / specla(j),
     >                      cbiol(i,j))
c
c update vegetation's physical characteristics
c
          plai(i,j)    = cbiol(i,j) * specla(j)
          biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
c
c update leaf and wood litter pools (we assume fire won't affect root pools)

          clitlm(i) = clitlm(i) * (1 - firefrac(i))
          clitls(i) = clitls(i) * (1 - firefrac(i))
          clitll(i) = clitll(i) * (1 - firefrac(i))
          clitwm(i) = clitwm(i) * (1 - firefrac(i))
          clitws(i) = clitws(i) * (1 - firefrac(i))
          clitwl(i) = clitwl(i) * (1 - firefrac(i))

 120    continue
c
c ---------------------------------------------------------------------
c * * * update annual npp, lai, and biomass * * *
c ---------------------------------------------------------------------
c
c adjust annual net ecosystem exchange (calculated in stats.f) 
c by loss of carbon to atmosphere due to biomass burning (fire)
c
        ayneetot(i) = ayneetot(i) - cdisturb(i)
c
c determine total ecosystem above-ground npp
c
        ayanpptot(i) = ayanpp(i,1)  + ayanpp(i,2) +
     >                 ayanpp(i,3)  + ayanpp(i,4) +
     >                 ayanpp(i,5)  + ayanpp(i,6) +
     >                 ayanpp(i,7)  + ayanpp(i,8) +
     >                 ayanpp(i,9)  + ayanpp(i,10) +
     >                 ayanpp(i,11) + ayanpp(i,12)
c
c update total canopy leaf area
c
        totlaiu(i) = plai(i,1)  + plai(i,2) +
     >               plai(i,3)  + plai(i,4) +
     >               plai(i,5)  + plai(i,6) +
     >               plai(i,7)  + plai(i,8)
c
        totlail(i) = plai(i,9)  + plai(i,10) +
     >               plai(i,11) + plai(i,12)
c
c update total biomass
c
        totbiou(i) = biomass(i,1) +
     >               biomass(i,2) +
     >               biomass(i,3) +
     >               biomass(i,4) +
     >               biomass(i,5) +
     >               biomass(i,6) +
     >               biomass(i,7) +
     >               biomass(i,8)
c
        totbiol(i) = biomass(i,9)  +
     >               biomass(i,10) +
     >               biomass(i,11) +
     >               biomass(i,12)
c
c ---------------------------------------------------------------------
c * * * update fractional cover and vegetation height parameters * * *
c ---------------------------------------------------------------------
c
c update fractional cover of forest and herbaceous canopies:
c 
        fu(i) = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
c
        fl(i) = totlail(i) / 1.0
c
c apply disturbances to fractional cover
c assume upper canopy burned 12%, lower canopy 100%

        fracupper = 0.12
        fraclower = 1.0

        fu(i) = fu(i) * (1 - firefrac(i)*fracupper)
        fl(i) = fl(i) * (1 - firefrac(i)*fraclower)

c
c constrain the fractional cover
c
        fu(i) = max (0.25, min (0.975, fu(i)))
        fl(i) = max (0.25, min (0.975, fl(i)))
c
c annual update upper canopy height parameters
c should be calculated based on vegetative fraction and not the
c average over the entire grid cell
c
        zbot(i,2) = 3.0
        ztop(i,2) = max(zbot(i,2) + 1.00, 2.50 *
     >                  totbiou(i) / fu(i) * 0.75)
c
c ---------------------------------------------------------------------
c * * * update stem area index and sapwood fraction * * *
c ---------------------------------------------------------------------
c
c estimate stem area index (sai) as a fraction of the lai
c
        sai(i,1) = 0.050 * totlail(i)
        sai(i,2) = 0.250 * totlaiu(i)
c
c estimate sapwood fraction of woody biomass
c
        sapspeed  = 25.0                        ! (m/day)
        trans     = 0.0025                      ! (2.5 mm/day) 
        saparea   = (trans / sapspeed)          ! m**2
c
        sapvolume = saparea * ztop(i,2) * 0.75  ! m**3
c
        denswood  = 400.0                       ! kg/m**3
c
        sapfrac(i) = min (0.50, max (0.05, sapvolume * denswood / wood))
c
       endif  ! check for crop existence 
c
 100  continue
c
c ---------------------------------------------------------------------
c * * * map out vegetation classes for this year * * *
c ---------------------------------------------------------------------
c
      call vegmap
c
c return to the main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine fire(iyear)
c ---------------------------------------------------------------------
c This subroutine written April 2010 by MMM and is based largely on globFIRM
c fire model used by LPJ. Refer to Thonicke 2004 for derivations and details.
c Previous fire routine saved in earlier IBIS versions, or see MMM.

      use comgrid
      use compar
      use comveg
      use comsoi
      use combcs,only: cmask

      implicit none

c
c local variables
c
      integer i, iyear
c
      real temp, s, daysinyear
c
      logical is_leap

c determine days in year - in case of leap year
c
      daysinyear = 365
 
      if (is_leap(iyear)) then
           daysinyear = 366
      endif
c
c begin global grid
c firefrac calculation taken from Thonicke 2004, globFIRM fire model. 
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle

         if (totalit(i).gt.fuelthresh) then 
               s           = firelength(i)/daysinyear 
               temp        = 0.45*(s-1)**3 + 2.83*(s-1)**2 + 2.96*(s-1) + 1.04
               firefrac(i) = s * exp((s-1)/temp)
         endif

 100  continue      

      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine vegmap
c ---------------------------------------------------------------------
c
      use comgrid
      use compar
      use comveg
      use combcs,only: cmask
c
      implicit none
c
c local variables
c
      integer 
     >     i, j,            ! loop indice
     >     domtree          ! dominant tree
c
      real maxlai,          ! maximum lai
     >     totlai,          ! total ecosystem lai
     >     grassfrac,       ! fraction of total lai in grasses
     >     treefrac,        ! fraction of total lai in trees
     >     treelai,         ! lai of trees
     >     shrublai,        ! lai of shrubs
     >     grasslai,        ! lai of grass
     >     ratio,
     >     cropbio
c
c classify vegetation cover into standard ibis vegetation classes 
c
c ---------------------------------------------------
c  1: tropical evergreen forest / woodland
c  2: tropical deciduous forest / woodland
c  3: temperate evergreen broadleaf forest / woodland
c  4: temperate evergreen conifer forest / woodland
c  5: temperate deciduous forest / woodland
c  6: boreal evergreen forest / woodland
c  7: boreal deciduous forest / woodland
c  8: mixed forest / woodland
c  9: savanna
c 10: grassland / steppe 
c 11: dense shrubland
c 12: open shrubland
c 13: tundra
c 14: desert 
c 15: polar desert / rock / ice
c 16: croplands
c ---------------------------------------------------
c
c begin global grid
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c determine total lai and tree, shrub, and grass fractions
c
        treelai   = totlaiu(i) 
        shrublai  = plai(i,9)  + plai(i,10)
        grasslai  = plai(i,11) + plai(i,12)
c
c crop biomass -- as used as an overriding condition for
c determining a vegetation class
c
        cropbio  = 0.
        do 105 j = scpft, ecpft
         cropbio   = cropbio + biomass(i,j)
 105    continue
c 
c
        totlai    = max (0.01, totlail(i) + totlaiu(i))
c
c determine dominant tree type by lai dominance
c
        domtree = 0
        maxlai = 0.0
c
        do 110 j = 1, 8
          if (plai(i,j).gt.maxlai) then
            domtree = j
            maxlai = plai(i,j)
          endif
 110    continue
c
c assign initial vegetation type
c
        vegtype0(i) = -999.99
c
c dominant type:  tropical broadleaf evergreen tree
c
        if (domtree.eq.1) then
          if (treelai.gt.2.5)         vegtype0(i) =  1.0  ! tropical evergreen forest / woodland
          if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  tropical broadleaf drought-deciduous tree
c
        if (domtree.eq.2) then
          if (treelai.gt.2.5)         vegtype0(i) =  2.0  ! tropical deciduous forest / woodland
          if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  warm-temperate broadleaf evergreen tree
c
        if (domtree.eq.3) then
          if (treelai.gt.2.5)         vegtype0(i) =  3.0  ! temperate evergreen broadleaf forest / woodland
          if (treelai.le.2.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  temperate conifer evergreen tree
c
        if (domtree.eq.4) then
          if (treelai.gt.1.5)         vegtype0(i) =  4.0  ! temperate evergreen conifer forest / woodland
          if (treelai.le.1.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  temperate broadleaf deciduous tree
c
        if (domtree.eq.5) then
          if (treelai.gt.1.5)         vegtype0(i) =  5.0  ! temperate deciduous forest / woodland
          if (treelai.le.1.5)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.0.5) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c dominant type:  boreal conifer evergreen tree
c
        if (domtree.eq.6)             vegtype0(i) =  6.0  ! boreal evergreen forest / woodland

       if (domtree.eq.6) then
         if (treelai.gt.1.0)         vegtype0(i) =  6.0  ! boreal evergreen forest / woodland
         if (treelai.le.1.0) then
           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
       endif
c
c dominant type:  boreal broadleaf cold-deciduous tree
c
        if (domtree.eq.7)             vegtype0(i) =  7.0  ! boreal deciduous forest / woodland

       if (domtree.eq.7) then
         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
         if (treelai.le.1.0) then
           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
       endif
c
c dominant type:  boreal conifer cold-deciduous tree
c
        if (domtree.eq.8)             vegtype0(i) =  7.0  ! boreal deciduous forest / woodland

       if (domtree.eq.8) then
         if (treelai.gt.1.0)         vegtype0(i) =  7.0  ! boreal deciduous forest / woodland
         if (treelai.le.1.0) then
           if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
           if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
         endif
       endif
c
c temperate/boreal forest mixtures
c
        if ((domtree.ge.4).and.(domtree.le.8)) then
          ratio = (plai(i,5) + plai(i,7) + plai(i,8)) / 
     >            (plai(i,4) + plai(i,5) + plai(i,6) + 
     >             plai(i,7) + plai(i,8))
          if (treelai.gt.1.0) then
            if ((ratio.gt.0.45).and.(ratio.lt.0.55)) vegtype0(i) = 8.
          endif
          if ((domtree.le.5).and.(treelai.le.1.0)) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c no tree is dominant
c
        if (domtree.eq.0) then
          if (treelai.gt.1.0)         vegtype0(i) =  9.0  ! savanna
          if (treelai.le.1.0) then
            if (grasslai.ge.shrublai) vegtype0(i) = 10.0  ! grassland
            if (shrublai.ge.grasslai) vegtype0(i) = 11.0  ! closed shrubland
          endif
        endif
c
c overriding vegtation classifications
c
        if (totlai.lt.1.0)            vegtype0(i) = 12.0  ! open shrubland
        if (totlai.le.0.4)            vegtype0(i) = 14.0  ! desert
c
c overriding climatic rules
c
        if (gdd5(i).lt.350.0) then
          if (totlai.ge.0.4)          vegtype0(i) = 13.0  ! tundra
          if (totlai.lt.0.4)          vegtype0(i) = 15.0  ! polar desert
        endif
c
        if (gdd0(i).lt.100.0)         vegtype0(i) = 15.0  ! polar desert
c
        if (cropbio .gt. 0.0)         vegtype0(i) = 16.0  ! croplands
c
 100  continue
c
c return to the main program
c
      return
      end
c
