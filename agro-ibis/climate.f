c last updated C. Kucharik 10.04.02
c
c  ####   #          #    #    #    ##     #####  ######
c #    #  #          #    ##  ##   #  #      #    #
c #       #          #    # ## #  #    #     #    #####
c #       #          #    #    #  ######     #    #
c #    #  #          #    #    #  #    #     #    #
c  ####   ######     #    #    #  #    #     #    ######
c
c ---------------------------------------------------------------------
      subroutine climanl
c ---------------------------------------------------------------------
c
c this subsroutine is only used to initialize growing degree days,
c coldest temp, and warmest temp at very beginning - provides a
c climate 'history' based on monthly mean values
c
c uses:
c
      use comgrid
      use compar
      use combcs
      use comveg
c
      implicit none
c
c Local variables
c
      integer it1w,      ! indice of previous month (interpolation)
     >        it2w,      ! indice of following month (interpolation)
     >      i,k,lda      ! loop indices
c
      real rwork,        ! work variable (1/ndaypm)
     >     dt,           ! used for interpolation
     >     dtemp,        ! interpolated temperature
     >     dtmin,        ! interpolated daily min temperature
     >     dtmax         ! interpolated daily max temperature
c
c Externals
c
      real gdd_natveg, gdd_crop  ! from gdd.f

c
c initialize values
c
      call const (gdd0,  npoi, 0.0)
      call const (gdd0c, npoi, 0.0)
      call const (gdd5,  npoi, 0.0)
      call const (gdd8,  npoi, 0.0)
      call const (gdd10, npoi, 0.0)
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c Find average annual precipitation in mm day-1 from monthly average
c Find average annual temperature in deg C from monthly average
c These are for phenology
c
        PPTavgann(i) = (xinprec(i,1) + xinprec(i,2) + xinprec(i,3) +
     >               xinprec(i,4) +  xinprec(i,5) +  xinprec(i,6) +
     >               xinprec(i,7) +  xinprec(i,8) +  xinprec(i,9) +
     >               xinprec(i,10) +  xinprec(i,11) + xinprec(i,12))/12.
c
        Tavgann(i) = (xint(i,1) + xint(i,2) + xint(i,3) +
     >               xint(i,4) +  xint(i,5) + xint(i,6) +
     >               xint(i,7) +  xint(i,8) + xint(i,9) +
     >               xint(i,10) + xint(i,11) + xint(i,12))/12.
c
c coldest monthly temperature (year 0) in deg c
c
        tc(i) = min (xint(i,1),  xint(i,2),  xint(i,3),
     >               xint(i,4),  xint(i,5),  xint(i,6),
     >               xint(i,7),  xint(i,8),  xint(i,9),
     >               xint(i,10), xint(i,11), xint(i,12))
c
c warmest monthly temperature (year 0) in deg c
c
        tw(i) = max (xint(i,1),  xint(i,2),  xint(i,3),
     >               xint(i,4),  xint(i,5),  xint(i,6),
     >               xint(i,7),  xint(i,8),  xint(i,9),
     >               xint(i,10), xint(i,11), xint(i,12))
c
        tcmin(i) = tc(i) + deltat(i)
c
 100  continue 
c
c interpolating climatological monthly input values to daily
c
      do 200 i = 1, npoi 
        if(cmask(i) == 0 ) cycle
c
        do 210 k = 1, 12
c
          rwork = 1. / float(ndaypm(k))
c
          do 220 lda = 1, ndaypm(k)
c
            if (float(lda).lt.float(ndaypm(k)+1)*0.5) then
              it1w = k - 1
              it2w = k
              dt   = (float(lda) - 0.5) * rwork + 0.5
            else
              it1w = k
              it2w = k + 1
              dt   = (float(lda) - 0.5) * rwork - 0.5
            end if
c
            if (it1w.lt. 1) it1w = 12
            if (it2w.gt.12) it2w = 1
c
            dtemp = xint(i,it1w) +
     >              dt * (xint(i,it2w) - xint(i,it1w))
            dtmin = xintmin(i,it1w) +
     >              dt * (xintmin(i,it2w) - xintmin(i,it1w))
            dtmax = xintmax(i,it1w) +
     >              dt * (xintmax(i,it2w) - xintmax(i,it1w))
c
c growing degree days, using deg c
c
            gdd0(i) = gdd0(i) + gdd_natveg(dtemp, 0.0)
            gdd5(i) = gdd5(i) + gdd_natveg(dtemp, 5.0)
c
c for crop phenological controls, impose an upper limit of
c of 30 C on gdd accumulation 
c
c
c         if (k .gt. 3 .and. k .lt. 10) then
            gdd0c(i) = gdd0c(i) + gdd_crop(dtemp, dtmax, dtmin, 0.0, 26.0)
            gdd8(i)  = gdd8(i)  + gdd_crop(dtemp, dtmax, dtmin, 8.0, 30.0)
            gdd10(i) = gdd10(i) + gdd_crop(dtemp, dtmax, dtmin, 10.0, 30.0)
c         endif
c
c
 220      continue
c
 210    continue
c
 200  continue
c
c call routine to determine pft existence arrays
c
      call existence

c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine climanl2
c ---------------------------------------------------------------------
c
c this subroutine updates the growing degree days, coldest temp, and
c warmest temp if monthly anomalies or daily values are used
c
c uses:
c
      use comgrid
      use compar
      use combcs
      use comveg
c
      implicit none
c
c local variables
c
      integer i            ! loop indice
c
      real zweigc,          ! 30-year e-folding time-avarage
     >     zweigw,          ! 30-year e-folding time-avarage
     >     rworkc,          ! 30-year e-folding time-avarage
     >     rworkw 
c
c calculate a 30-year e-folding time-avarage
c
      zweigc = exp(-1./30.)
      zweigw = exp(-1./30.)
c
      rworkc = 1. - zweigc
      rworkw = 1. - zweigw
c
c update critical climatic parameters with running average
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
        tc(i) = zweigc * tc(i) + rworkc * tcthis(i)
        tw(i) = zweigw * tw(i) + rworkw * twthis(i)
c
        tcmin(i) = tc(i) + deltat(i)
c
        gdd0(i) = zweigc * gdd0(i) +
     >            rworkc * gdd0this(i)
c
        gdd5(i) = zweigc * gdd5(i) +
     >            rworkc * gdd5this(i)
c
        gdd0c(i) = zweigc * gdd0c(i) +
     >            rworkc * gdd0cthis(i)
c
        gdd8(i) = zweigc * gdd8(i) +
     >            rworkc * gdd8this(i)
c
        gdd10(i) = zweigc * gdd10(i) +
     >            rworkc * gdd10this(i)
c
 100  continue
c
      call existence
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine existence
c ---------------------------------------------------------------------
c
c this routine determines which plant functional types (pft's) are allowed
c to exist in each gridcell, based on a simple set of climatic criteria
c
c the logic here is based on the biome3 model of haxeltine and prentice
c
c plant functional types:
c
c 1)  tropical broadleaf evergreen trees
c 2)  tropical broadleaf drought-deciduous trees
c 3)  warm-temperate broadleaf evergreen trees
c 4)  temperate conifer evergreen trees
c 5)  temperate broadleaf cold-deciduous trees
c 6)  boreal conifer evergreen trees
c 7)  boreal broadleaf cold-deciduous trees
c 8)  boreal conifer cold-deciduous trees
c 9)  evergreen shrubs
c 10) deciduous shrubs
c 11) warm (c4) grasses
c 12) cool (c3) grasses
c 13) c3 crop (soybean)
c 14) c4 crop (corn)
c 15) c3 crop (wheat)
c
      use comgrid
      use compar
      use comveg
      use comsum
      use combcs
      use comcrop
      use comnitr
      use compft
      use comwork
c
      implicit none
c
c Local variables
c
      integer i,j,        ! loop index
     >        inveg        
c
c ---------------------------------------------------------------------
c
      do 100 i = 1, npoi
        if(cmask(i) == 0 ) cycle
c
c determine which plant types can exist in a given gridcell
c
        exist(i,1)  = 0.
        exist(i,2)  = 0.
        exist(i,3)  = 0.
        exist(i,4)  = 0.
        exist(i,5)  = 0.
        exist(i,6)  = 0.
        exist(i,7)  = 0.
        exist(i,8)  = 0.
        exist(i,9)  = 0.
        exist(i,10) = 0.
        exist(i,11) = 0.
        exist(i,12) = 0.
        exist(i,13) = 0.
        exist(i,14) = 0.
        exist(i,15) = 0.
c
c 1) tropical broadleaf evergreen trees
c
c  - tcmin > 0.0
c
*       if (tcmin(i).gt.0.0)           exist(i,1) = 1.0
c
c 2) tropical broadleaf drought-deciduous trees
c
c  - tcmin > 0.0
c
*        if (tcmin(i).gt.0.0)           exist(i,2) = 1.0
c
c 3) warm-temperate broadleaf evergreen trees
c
c  - tcmin <   0.0 and
c  - tcmin > -10.0
c
*        if ((tcmin(i).lt.0.0).and.
*     >      (tcmin(i).gt.-10.0))       exist(i,3) = 1.0
c
c 4) temperate conifer evergreen trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
*        if ((tcmin(i).lt.0.0).and.
*     >      (tcmin(i).gt.-45.0).and.
*     >      (gdd5(i).gt.1200.0))       exist(i,4) = 1.0
c
c 5) temperate broadleaf cold-deciduous trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
*        if ((tcmin(i).lt.0.0).and.
*     >      (tcmin(i).gt.-45.0).and.
*     >      (gdd5(i).gt.1200.0))       exist(i,5) = 1.0
c
c 6) boreal conifer evergreen trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
*       if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
*     >       (tcmin(i).gt.-57.5).and.
*     >       (gdd5(i).gt.350.0))       exist(i,6) = 1.0
c
c 7) boreal broadleaf cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
*       if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
*     >       (tcmin(i).gt.-57.5).and.
*     >       (gdd5(i).gt.350.0))       exist(i,7) = 1.0
c
c 8) boreal conifer cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - gdd5  >  350.0
c
*        if (((tcmin(i).lt.-45.0).or.(gdd5(i).lt.1200.0)).and.
*     >       (gdd5(i).gt.350.0))       exist(i,8) = 1.0
c
c 9) evergreen shrubs
c
c  - gdd0 > 100.0
c
*        if (gdd0(i).gt.100.0)          exist(i,9) = 1.0
c
c 10) deciduous shrubs
c
c  - gdd0 > 100.0
c
*        if (gdd0(i).gt.100.0)          exist(i,10) = 1.0
c
c 11) warm (c4) grasses
c
c  - tw   >  22.0 and
c  - gdd0 > 100.0
c
*        if ((tw(i).gt.22.0).and.
*     >      (gdd0(i).gt.100.0))        exist(i,11) = 1.0
c
c 12) cool (c3) grasses
c
c  - gdd0 > 100.0
c
*        if (gdd0(i).gt.100.0)          exist(i,12) = 1.0
        
**** DTP 2001/06/07: Modified version of above code reads in PFT
*    existence criteria from external parameter file "params.veg"
*    These are copied here for reference.... 
*------------------------------------------------------------------
*  TminL    TminU    Twarm    GDD    PFT
*------------------------------------------------------------------
*    0.0   9999.0   9999.0   9999  !   1
*    0.0   9999.0   9999.0   9999  !   2
*  -10.0      0.0   9999.0   9999  !   3
*  -45.0      0.0   9999.0   1200  !   4
*  -45.0      0.0   9999.0   1200  !   5
*  -57.5    -45.0   9999.0    350  !   6
*  -57.5    -45.0   9999.0    350  !   7
* 9999.0    -45.0   9999.0    350  !   8
* 9999.0   9999.0   9999.0    100  !   9
* 9999.0   9999.0   9999.0    100  !  10
* 9999.0   9999.0     22.0    100  !  11
* 9999.0   9999.0   9999.0    100  !  12
*------------------------------------------------------------------

c 1) tropical broadleaf evergreen trees
c
c  - tcmin > 0.0
c
        if (tcmin(i).gt.TminL(1))      exist(i,1) = 1.0
c
c 2) tropical broadleaf drought-deciduous trees
c
c  - tcmin > 0.0
c
        if (tcmin(i).gt.TminL(2))      exist(i,2) = 1.0
c
c 3) warm-temperate broadleaf evergreen trees
c
c  - tcmin <   0.0 and
c  - tcmin > -10.0
c
        if ((tcmin(i).lt.TminU(3)).and.
     >      (tcmin(i).gt.TminL(3)))    exist(i,3) = 1.0
c
c 4) temperate conifer evergreen trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
        if ((tcmin(i).lt.TminU(4)).and.
     >      (tcmin(i).gt.TminL(4)).and.
     >      (gdd5(i).gt.GDD(4)))       exist(i,4) = 1.0
c
c 5) temperate broadleaf cold-deciduous trees
c
c  - tcmin <    0.0 and
c  - tcmin >  -45.0 and
c  - gdd5  > 1200.0
c
        if ((tcmin(i).lt.TminU(5)).and.
     >      (tcmin(i).gt.TminL(5)).and.
     >      (gdd5(i).gt.GDD(5)))       exist(i,5) = 1.0
c
c 6) boreal conifer evergreen trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
        if (((tcmin(i).lt.TminU(6)).or.
     >      (gdd5(i).lt.GDD(4))).and.
     >      (tcmin(i).gt.TminL(6)).and.
     >      (gdd5(i).gt.GDD(6)))       exist(i,6) = 1.0
c
c 7) boreal broadleaf cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - tcmin >  -57.5 and
c  - gdd5  >  350.0
c
        if (((tcmin(i).lt.TminU(7)).or.
     >      (gdd5(i).lt.GDD(5))).and.
     >      (tcmin(i).gt.TminL(7)).and.
     >      (gdd5(i).gt.GDD(7)))       exist(i,7) = 1.0
c
c 8) boreal conifer cold-deciduous trees
c
c  - tcmin <  -45.0 or gdd5 < 1200.0, and
c  - gdd5  >  350.0
c
        if (((tcmin(i).lt.TminU(8)).or.
     >      (gdd5(i).lt.TminL(4))).and.
     >      (gdd5(i).gt.GDD(8)))       exist(i,8) = 1.0
c
c 9) evergreen shrubs
c
c  - gdd0 > 100.0
c
        if (gdd0(i).gt.GDD(9))         exist(i,9) = 1.0
c
c 10) deciduous shrubs
c
c  - gdd0 > 100.0
c
        if (gdd0(i).gt.GDD(10))        exist(i,10) = 1.0
c
c 11) warm (c4) grasses
c
c  - tw   >  22.0 and
c  - gdd0 > 100.0
c
        if ((tw(i).gt.Twarm(11)).and.
     >      (gdd0(i).gt.GDD(11)))      exist(i,11) = 1.0
c
c 12) cool (c3) grasses
c
c  - gdd0 > 100.0
c
        if (gdd0(i).gt.GDD(12))        exist(i,12) = 1.0

c
c
c == C. Kucharik 6.12.01 ==
c
c if override natural vegetation competition (overveg = 1)
c this code is used to override existence parameterizations for potential
c vegetation distribution based on climatic constraints. Instead
c we only allow PFTs to compete in each grid cell
c based on land cover dataset and classification found in that region 
c override those pfts that are not desired but might have exist = 1.0
c from above initialization - this essentially limits vegetation competition
c during spin-up periods so that vegetation growing there is confined to
c what is typically observed today (potential vegetation).  If doing 
c climate change scenarios, overveg should be set to 0 so full 
c vegetation dynamics are used, if desired. 
c
        if (overveg .eq. 1) then
c
           inveg = nint(xinveg(i))
c          inveg = 6 
c          
c wherever vegetation is allowed to exist, check the rules from above 
c
c tropical deciduous
c
          if (inveg .eq. 2) then
            exist(i,2) = 1.0
c
            exist(i,1) = 0.0
            do 70 j = 3,15
              exist(i,j) = 0.0
 70         continue

c
c temperate conifers
c
          else if (inveg .eq. 4) then
            exist(i,4) = 1.0
c
            do 71 j = 1,3
              exist(i,j) = 0.0
 71         continue
c
            do 72 j = 5,15
              exist(i,j) = 0.0
 72         continue
c
c temperate deciduous
c
          else if (inveg .eq. 5) then
            exist(i,5) = 1.0
c
            do 73 j = 1,4
              exist(i,j) = 0.0
 73         continue
c
            do 74 j = 6,15
              exist(i,j) = 0.0
 74         continue
c
c boreal conifers
c
          else if (inveg .eq. 6) then
            exist(i,6) = 1.0
c
            do 75 j = 1,5
              exist(i,j) = 0.0
 75         continue
c
            do 76 j = 7,15
              exist(i,j) = 0.0
 76         continue
c
c boreal deciduous 
c
          else if (inveg .eq. 7) then
            do 77 j = 7, 8
              exist(i,j) = 1.0
 77         continue
c
            do 78 j = 1,6
              exist(i,j) = 0.0
 78         continue
c
            do 79 j = 9,15
              exist(i,j) = 0.0
 79         continue
c
c mixed forest exception:
c
c let existence rules determine whether temperate or boreal species
c can exist at this location
c 
          else if (inveg .eq. 8) then
c            do 69 j = 4, 8
c              exist(i,j) = 1.0
c 69         continue
c
            do 81 j = 1,3
              exist(i,j) = 0.0
 81         continue
c
            do 82 j = 9,15
              exist(i,j) = 0.0
 82         continue
c
c savanna
c
          else if (inveg .eq. 9) then
              exist(i,5) = 1.0
              exist(i,11) = 1.0
              exist(i,12) = 1.0
c
            do 83 j = 1, 4
              exist(i,j) = 0.0
 83         continue
c
            do 84 j = 6, 10 
              exist(i,j) = 0.0
 84         continue
c
            do 85 j = 13,15
              exist(i,j) = 0.0
 85         continue
c
c grassland
c
          else if (inveg .eq. 10) then
              exist(i,11) = 1.0
              exist(i,12) = 1.0
c
            do 86 j = 1, 10 
              exist(i,j) = 0.0
 86         continue
c
            do 87 j = 13, 15 
              exist(i,j) = 0.0
 87         continue
c
c dense shrubland 
c
          else if (inveg .eq. 11) then
              exist(i,9)  = 1.0
              exist(i,10) = 1.0
              exist(i,11) = 1.0
              exist(i,12) = 1.0
c
            do 88 j = 1, 8 
              exist(i,j) = 0.0
 88         continue
c
            do 89 j = 13, 15 
              exist(i,j) = 0.0
 89         continue
c
c open shrubland
c
          else if (inveg .eq. 12) then
c
            do 90 j = 9, 12 
              exist(i,j) = 1.0
 90         continue
c
            do 91 j = 1, 8 
              exist(i,j) = 0.0
 91         continue
c
            do 92 j = 13, 15 
              exist(i,j) = 0.0
 92         continue
c
c tundra 
c
          else if (inveg .eq. 13) then
c
            do 93 j = 9, 12 
              exist(i,j) = 1.0
 93         continue
c
            do 94 j = 1, 8 
              exist(i,j) = 0.0
 94         continue
c
            do 95 j = 13, 15 
              exist(i,j) = 0.0
 95         continue
          endif
c
         endif ! override original natural vegetation competition
c
c == C. Kucharik ==
c cropping systems - are currently planted everywhere
c
          icropsum(i) = isoybean + imaize + iwheat + irotation
c
          if (icropsum(i) .gt. 0.0) then  
            do 80 j = 1, 12
              exist(i,j) = 0.0
 80         continue
          endif
c
c 13) c3 crop - soybean
c
c if crops exist at a grid cell, natural vegetation
c is not allowed to exist  
c
          if (isoybean .eq. 1)          exist(i,13) = 1.0 
c
c 14) c4 crop - corn
c
          if (imaize .eq. 1)            exist(i,14) = 1.0
c
c 15) c3 crop - wheat (spring and winter varieties) 
c
          if (iwheat .gt. 0)            exist(i,15) = 1.0
c
c
 100  continue
c
      return
      end
