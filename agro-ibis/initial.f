c last update 1.26.05 C. Kucharik for Green-Ampt
c
c    #    #    #     #     #####     #      ##    #
c    #    ##   #     #       #       #     #  #   #
c    #    # #  #     #       #       #    #    #  #
c    #    #  # #     #       #       #    ######  #
c    #    #   ##     #       #       #    #    #  #
c    #    #    #     #       #       #    #    #  ######
c
c ---------------------------------------------------------------------
      subroutine initial (isimveg, irestart, iyrlast)
c ---------------------------------------------------------------------
c
      use combcs,only: cmask

      implicit none    
c
c Arguments (input)     
c
      integer isimveg,          ! 0 = static veg, 1 = dynamic veg, 
     >                          ! 2 = dynamic veg with cold start
     >        irestart,         ! 0 = initial run, 1 = restart run
     >     iyrlast              ! last year of previous run (for restart)
c
      if (irestart .eq. 0) then
        call coldstart
      else
        call restart (iyrlast)
      end if

c
c initialize physical consts, dimensions, unit numbers, lsx model
c
      call inisurf(irestart)
c
c initialize snow model
c
      call inisnow
c
c initialize soil model
c
      call inisoil (irestart)
c
c initialize vegetation parameters
c
      call iniveg (isimveg, irestart)
c
c initialize variables for time averaging
c
      call inisum
c
c return to main program
c
      return
      end

c
c ---------------------------------------------------------------------
      subroutine inisurf(irestart)
c ---------------------------------------------------------------------
c
c does initialization for model
c
      use comgrid
      use compar
      use comatm
      use comhyd
      use comsum
      use comveg
c
      implicit none
c
c Arguments (input)     
c
      integer irestart          ! 0 = initial run, 1 = restart run
c
c local variables
c
      integer i,
     >        j       ! loop indice
c
c set physical constants (mks)
c
      stef  = 5.67e-8 
      vonk  = 0.4
      grav  = 9.80616
      tmelt = 273.16
      hvap  = 2.5104e+6
      hfus  = 0.3336e+6
      hsub  = hvap + hfus
      ch2o  = 4.218e+3
      cice  = 2.106e+3
      cair  = 1.004640000e+3
      cvap  = 1.810000000e+3
      rair  = 287.04
      rvap  = 461.0
      cappa = rair / cair
      rhow  = 1.0e+3
c
      call const (vzero, npoi, 0.0)
c
c specify the epsilon value for the model
c
      epsilon = 1.0e-7
c
c initialize integer variables (can't use const for this)
c
c wet day / dry day flag initialized to dry day (0)
c
      do 100 i = 1, npoi 
        iwet(i) = 0
        do 105 j = 1,31
          iwetday(i,j) = 0
          precipday(i,j) = 0
 105    continue 
 100  continue
c
c zero flux arrays, and global diagnostic arrays
c
      call const (asurd, npoi*nband, 0.0)
      call const (asuri, npoi*nband, 0.0)
c
      call const (totcondub, npoi, 0.0)
      call const (totconduc, npoi, 0.0)
      call const (totcondls, npoi, 0.0)
      call const (totcondl3, npoi, 0.0)
      call const (totcondl4, npoi, 0.0)
      call const (totcondc3, npoi, 0.0)
      call const (totcondc4, npoi, 0.0)
c
      call const (ginvap, npoi, 0.0)
      call const (gsuvap, npoi, 0.0)
      call const (gtrans, npoi, 0.0)
      call const (grunof, npoi, 0.0)
      call const (gdrain, npoi, 0.0)
      call const (gdrain_layer, npoi*nsoilay, 0.0)
c
c initialize vegetation prognostic variables
c
c initialize all temperature fields to 10 degrees C
c
      call const (tu,    npoi, 283.16)
      call const (ts,    npoi, 283.16)
      call const (tl,    npoi, 283.16)
c
c initialize weather generator 'memory'
c
      call const (xstore, npoi*3, 0.0)
c
c initialize temperature of lower canopy buried by
c snow to 0 degrees C
c
      call const (tlsub, npoi, 273.16)
c
c initialize canopy air conditions (used in turvap)
c
      call const (t12, npoi, 283.16)
      call const (t34, npoi, 283.16)
c
      call const (q12, npoi, 0.0)
      call const (q34, npoi, 0.0)
c
c initialize all co2 concentrations (mol/mol)
c
      call const (ciub, npoi, 350.0e-06) 
      call const (ciuc, npoi, 350.0e-06) 
      call const (cils, npoi, 350.0e-06) 
      call const (cil3, npoi, 350.0e-06) 
      call const (cil4, npoi, 350.0e-06) 
      call const (cic3, npoi, 350.0e-06)
      call const (cic4, npoi, 350.0e-06)
c
      call const (csub, npoi, 350.0e-06) 
      call const (csuc, npoi, 350.0e-06) 
      call const (csls, npoi, 350.0e-06) 
      call const (csl3, npoi, 350.0e-06) 
      call const (csl4, npoi, 350.0e-06) 
      call const (csc3, npoi, 350.0e-06)
      call const (csc4, npoi, 350.0e-06)
c
c initialize stomatal conductance (mol-h2o/m**2/sec)
c
      call const (gsub, npoi, 0.5) 
      call const (gsuc, npoi, 0.5) 
      call const (gsls, npoi, 0.5) 
      call const (gsl3, npoi, 0.5) 
      call const (gsl4, npoi, 0.5) 
      call const (gsc3, npoi, 0.5)
      call const (gsc4, npoi, 0.5)
c
c initialize soil biogeochemistry variables
c
      if (irestart .eq. 0) then
         call const (clitlm,   npoi, 0.0) 
         call const (clitls,   npoi, 0.0) 
         call const (clitll,   npoi, 0.0) 
         call const (clitrm,   npoi, 0.0) 
         call const (clitrs,   npoi, 0.0) 
         call const (clitrl,   npoi, 0.0) 
         call const (clitwm,   npoi, 0.0) 
         call const (clitws,   npoi, 0.0) 
         call const (clitwl,   npoi, 0.0) 
         call const (totcmic,  npoi, 0.0) 
         call const (csoislop, npoi, 0.0) 
         call const (csoislon, npoi, 0.0) 
         call const (csoipas,  npoi, 0.0) 
      end if
      call const (totlit,   npoi, 0.0) !
      call const (totnlit,  npoi, 0.0) !
      call const (totfall,  npoi, 0.0) !
      call const (totalit,  npoi, 0.0) !
      call const (totrlit,  npoi, 0.0) !
      call const (totanlit, npoi, 0.0) !
      call const (totrnlit, npoi, 0.0) !
      call const (totcsoi,  npoi, 0.0) !
      call const (totnmic,  npoi, 0.0) !
      call const (tco2mic,  npoi, 0.0) !
      call const (tnpptot,  npoi, 0.0) !
      call const (tneetot,  npoi, 0.0) !
      call const (tnmin,    npoi, 0.0) !
c
c initialize carbon lost to atmosphere due
c to biomass burning
c
      call const (cdisturb, npoi, 0.0)
c
c initialize phenology flags
c
      if (irestart .eq. 0) then
         call const (tempu,  npoi, 1.0)
         call const (templs,  npoi, 1.0)
         call const (greenfracl3, npoi, 1.0)
         call const (greenfracl4, npoi, 1.0)
c
         call logicf (onflagu, npoi)
         call logicf (onflagl3, npoi)
         call logicf (onflagl4, npoi)
         call logict (offflagu, npoi)
         call logict (offflagl3, npoi)
         call logict (offflagl4, npoi)
      end if
c
c initialize water and snow interception fractions
c
      call const (wliqu, npoi, 0.0)
      call const (wliqs, npoi, 0.0)
      call const (wliql, npoi, 0.0)
c
      call const (wsnou, npoi, 0.0)
      call const (wsnos, npoi, 0.0)
      call const (wsnol, npoi, 0.0)
c
      call const (su, npoi, 0.0)
      call const (ss, npoi, 0.0)
      call const (sl, npoi, 0.0)
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inisum
c ---------------------------------------------------------------------
c CD
c does initialization for time averaging
c
      use comgrid
      use compar
      use comhyd
      use comsno
      use comsoi
      use comsum
      use comveg
c
      implicit none
c
c local variables
c
      integer i, k        ! loop indices
c
c initialize total water content in soil+snow+vegetation (for mass conservation check)
c
      do 20 i = 1, npoi
        wtot(i) = (wliqu(i)+wsnou(i)) * fu(i) * 2.0 * lai(i,2) +
     >            (wliqs(i)+wsnos(i)) * fu(i) * 2.0 * sai(i,2) +
     >            (wliql(i)+wsnol(i)) * fl(i) * 2.0 *
     >            (lai(i,1) + sai(i,1)) * (1. - fi(i))
c
        wtot(i) = wtot(i) + wpud(i) + wipud(i)
c
        do 10 k = 1, nsoilay
          wtot(i) = wtot(i) +
     >              poros(i,k)*wsoi(i,k)*(1.-wisoi(i,k))*hsoi(k)*rhow +
     >              poros(i,k)*wisoi(i,k)*hsoi(k)*rhow
 10     continue
c
        do 20 k = 1, nsnolay
          wtot(i) = wtot(i) + fi(i)*rhos*hsno(i,k)
 20     continue
c
c Daily means
c
      call const (adrain,     npoi, 0.0) 
      call const (adsnow,     npoi, 0.0) 
      call const (adaet,      npoi, 0.0) 
      call const (adtrans,    npoi, 0.0) 
      call const (adevap,     npoi, 0.0) 
      call const (adtratio,   npoi, 0.0) 
      call const (adtrunoff,  npoi, 0.0) 
      call const (adsrunoff,  npoi, 0.0) 
      call const (addrainage, npoi, 0.0) 
      call const (addrainage_layer, npoi*nsoilay, 0.0) 
      call const (adrh,       npoi, 0.0) 
      call const (adud,       npoi, 0.0) 
      call const (adsnod,     npoi, 0.0) 
      call const (adsnof,     npoi, 0.0) 
      call const (adwsoi,     npoi, 0.0) 
      call const (adtsoi,     npoi, 0.0) 
      call const (adwisoi,    npoi, 0.0)
      call const (adtlaysoi,  npoi, 0.0) 
      call const (adwlaysoi,  npoi, 0.0)
      call const (adwsoic,    npoi, 0.0) 
      call const (adtsoic,    npoi, 0.0) 
      call const (adco2mic,   npoi, 0.0) 
      call const (adco2root,  npoi, 0.0) 
      call const (decompl,    npoi, 0.0) 
      call const (decomps,    npoi, 0.0) 
      call const (adnmintot,  npoi, 0.0) 
c
c monthly mean quanties
c
      call const (amtemp,     npoi, 0.0)
      call const (amrain,     npoi, 0.0) 
      call const (amsnow,     npoi, 0.0) 
      call const (amaet,      npoi, 0.0) 
      call const (amtrunoff,  npoi, 0.0) 
      call const (amsrunoff,  npoi, 0.0) 
      call const (amdrainage, npoi, 0.0) 
      call const (amdrainage_layer, npoi*nsoilay, 0.0) 
      call const (amtotnleach,npoi*nsoilay, 0.0) 
      call const (amno3leach, npoi, 0.0) 
      call const (amcloud,    npoi, 0.0)
      call const (amqa,       npoi, 0.0)
      call const (amrh,       npoi, 0.0) 
      call const (amsolar,    npoi, 0.0) 
      call const (amsolarnet, npoi, 0.0) 
      call const (amirup,     npoi, 0.0) 
      call const (amirdown,   npoi, 0.0) 
      call const (amsens,     npoi, 0.0) 
      call const (amlatent,   npoi, 0.0) 
      call const (amnetrad,   npoi, 0.0) 
      call const (amlaiu,     npoi, 0.0) 
      call const (amlail,     npoi, 0.0) 
      call const (amglail,    npoi, 0.0) 
      call const (amblail,    npoi, 0.0) 
      call const (amtsoi,     npoi, 0.0) 
      call const (amwsoi,     npoi, 0.0) 
      call const (amwisoi,    npoi, 0.0)
      call const (amvwc,      npoi, 0.0)
      call const (amawc,      npoi, 0.0)
      call const (amsnod,     npoi, 0.0) 
      call const (amsnof,     npoi, 0.0)
      call const (amco2mic,   npoi, 0.0) 
      call const (amco2root,  npoi, 0.0) 
      call const (amnmintot,  npoi, 0.0) 
c
      call const (amnpp,      npoi*npft, 0.0)     
c
c Annual mean quantities
c
      call const (aysolar,    npoi, 0.0) 
      call const (aysolarnet, npoi, 0.0) 
      call const (ayirup,     npoi, 0.0) 
      call const (ayirdown,   npoi, 0.0) 
      call const (aysens,     npoi, 0.0) 
      call const (aylatent,   npoi, 0.0) 
      call const (aynetrad,   npoi, 0.0) 
      call const (ayprcp,     npoi, 0.0) 
      call const (ayaet,      npoi, 0.0) 
      call const (aytrans,    npoi, 0.0) 
      call const (aytrunoff,  npoi, 0.0) 
      call const (aysrunoff,  npoi, 0.0) 
      call const (aydrainage, npoi, 0.0) 
      call const (aywsoi,     npoi, 0.0) 
      call const (aywisoi,    npoi, 0.0)
      call const (aytsoi,     npoi, 0.0) 
      call const (ayvwc,      npoi, 0.0)
      call const (ayawc,      npoi, 0.0)
      call const (aystresstu, npoi, 0.0)
      call const (aystresstl, npoi, 0.0)
      call const (firefac,    npoi, 0.0)
      call const (firefac,    npoi, 0.0)
      call const (ayco2mic,   npoi, 0.0) 
      call const (ayco2root,  npoi, 0.0) 
      call const (ayrootbio,  npoi, 0.0) 
      call const (aynmintot,  npoi, 0.0) 
      call const (ayalit,     npoi, 0.0)
      call const (ayblit,     npoi, 0.0)
      call const (aycsoi,     npoi, 0.0)
      call const (aycmic,     npoi, 0.0)
      call const (ayanlit,    npoi, 0.0)
      call const (aybnlit,    npoi, 0.0)
      call const (aynsoi,     npoi, 0.0)
c
      call const (aygpp,      npoi*npft, 0.0)     
      call const (aygpp,      npoi*npft, 0.0)    
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inisnow
c ---------------------------------------------------------------------
c
c does initialization for snow model
c
      use comgrid
      use compar
      use comsno
c
      implicit none
c
c rhos is density of snow
c
      rhos = 0.15 * rhow
c
c consno is thermal conductivity of snow
c
      consno = 0.20
c
c hsnotop is "adaptive-grid" thickness of top snow layer
c
      hsnotop = 0.1
c
c hsnomin is minimum total snow thickness. total thickness
c is constrained to hsnomin for less than 100% cover. (hsnomin
c should be ge nsnolay*hsnotop for vadapt to work properly.)
c
      hsnomin = max (0.3, nsnolay * hsnotop)
c
c fimin and fimax are minimum and maximum snowcover fractions
c
      fimin = 0.00002 * (dtime / 1800.) * (0.3 / hsnomin)
      fimax = 1.000
c
c z0sno is roughness lenth of snow cover
c
      z0sno = 0.0005
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inisoil (irestart)
c ---------------------------------------------------------------------
c
c does initialization for soil database
c
      use comgrid
      use compar
      use comsoi
      use comsno
      use comsum
      use comatm
      use comcrop
      use comnitr
      use comtex
      use comwork, only: cdummy
c
      implicit none
c
c Arguments (input)     
c
      integer irestart          ! 0 = initial run, 1 = restart run
c
c local variables
c
      integer i,
     >        k,
     >        l,
*    >        ndat,           ! number of textural classes
     >        msand,          ! % of sand in grid point
     >        mclay,          ! % of clay in grid point
     >        lmin,           ! closest textural class from texture of 
     >                        ! grid point
     >        textcls,         ! textural class assignment (1..11)
     >        domtexti
c
      real fsand,             ! fraction of sand in grid point
     >     fsilt,             ! fraction of silt in grid point
     >     fclay,             ! fraction of clay in grid point
     >     forganic           ! fraction of organic matter in soil

*      parameter (ndat=11)
c
      real 
*     >     texdat(3,ndat),    ! % of sand, silt, clay in each textural class
*     >     porosdat(ndat),    ! porosity in fraction of soil depth
*     >     sfielddat(ndat),   ! field capacity in fraction of soil depth
*     >     swiltdat(ndat),    ! wilting point in fraction of soil depth
*     >     bexdat(ndat),      ! 'b' exponent in moisture release equation
*     >     suctiondat(ndat),  ! saturated (air entry) potential (m-h2o)
*     >     hydrauldat(ndat),  ! saturated hydraulic conductivity (m s-1)
     >     xdat(ndat),        ! % of sand in each textural class
     >     ydat(ndat),        ! % of silt in each textural class
     >     zdat(ndat),		  ! % of clay in each textural class
     >     cpwfdat(ndat)      ! Green-Ampt, capillary pressure at wetting front (m)
c
c
c Rawls et al. (1992) soil properties data
c
c      ------------------
c       sand  silt  clay
c      ------------------
c
*      data texdat /
*     >  0.92, 0.05, 0.03,  ! sand
*     >  0.81, 0.12, 0.07,  ! loamy sand
*     >  0.65, 0.25, 0.10,  ! sandy loam
*     >  0.42, 0.40, 0.18,  ! loam
*     >  0.20, 0.65, 0.15,  ! silt loam
*     >  0.60, 0.13, 0.27,  ! sandy clay loam
*     >  0.32, 0.34, 0.34,  ! clay loam
*     >  0.09, 0.58, 0.33,  ! silty clay loam
*     >  0.53, 0.07, 0.40,  ! sandy clay
*     >  0.10, 0.45, 0.45,  ! silty clay
*     >  0.20, 0.20, 0.60   ! clay
*     >  /
c
c porosity (fraction)
c
*      data porosdat /
*     >  0.437,             ! sand
*     >  0.437,             ! loamy sand
*     >  0.453,             ! sandy loam
*     >  0.463,             ! loam
*     >  0.501,             ! silt loam
*     >  0.398,             ! sandy clay loam
*     >  0.464,             ! clay loam
*     >  0.471,             ! silty clay loam
*     >  0.430,             ! sandy clay
*     >  0.479,             ! silty clay
*     >  0.475              ! clay
*     >  /
c
c field capacity (fraction)
c
*      data sfielddat /
*     >  0.091,             ! sand
*     >  0.125,             ! loamy sand
*     >  0.207,             ! sandy loam
*     >  0.270,             ! loam
*     >  0.330,             ! silt loam
*     >  0.255,             ! sandy clay loam
*     >  0.318,             ! clay loam
*     >  0.366,             ! silty clay loam
*     >  0.339,             ! sandy clay
*     >  0.387,             ! silty clay
*     >  0.396              ! clay
*     >  /
c
c wilting point (fraction)
c
*      data swiltdat /
*     >  0.033,             ! sand
*     >  0.055,             ! loamy sand
*     >  0.095,             ! sandy loam
*     >  0.117,             ! loam
*     >  0.133,             ! silt loam
*     >  0.148,             ! sandy clay loam
*     >  0.197,             ! clay loam
*     >  0.208,             ! silty clay loam
*     >  0.239,             ! sandy clay
*     >  0.250,             ! silty clay
*     >  0.272              ! clay
*     >  /
c
c "b" exponent for the Campbell moisture-release equation
c
*      data bexdat /
*     >  1.7,               ! sand
*     >  2.1,               ! loamy sand
*     >  3.1,               ! sandy loam
*     >  4.5,               ! loam
*     >  4.7,               ! silt loam
*     >  4.0,               ! sandy clay loam
*     >  5.2,               ! clay loam
*     >  6.6,               ! silty clay loam
*     >  6.0,               ! sandy clay
*     >  7.9,               ! silty clay
*     >  7.6                ! clay
*     >  /
c
c saturated (air entry) potential (m-h2o)
c
*      data suctiondat /
*     >  0.070,             ! sand
*     >  0.090,             ! loamy sand
*     >  0.150,             ! sandy loam
*     >  0.110,             ! loam
*     >  0.210,             ! silt loam
*     >  0.280,             ! sandy clay loam
*     >  0.260,             ! clay loam
*     >  0.330,             ! silty clay loam
*     >  0.290,             ! sandy clay
*     >  0.340,             ! silty clay
*     >  0.370              ! clay
*     >  /
c
c saturated hydraulic conductivity (m s-1)
c
*      data hydrauldat /
*     >  5.8330e-05,        ! sand
*     >  1.6972e-05,        ! loamy sand
*     >  7.1944e-06,        ! sandy loam
*     >  3.6667e-06,        ! loam
*     >  1.8889e-06,        ! silt loam
*     >  1.1944e-06,        ! sandy clay loam
*     >  6.3889e-07,        ! clay loam
*     >  4.1667e-07,        ! silty clay loam
*     >  3.3333e-07,        ! sandy clay
*     >  2.5000e-07,        ! silty clay
*     >  1.6667e-07         ! clay
*     >  /
c
c
c --------------------------------------------------------------------------
c added for Green-Ampt infiltration model
c capillary pressure at the wetting front for Green-Ampt infiltration model
c units are meters of head (from Rawls et al. 1983, ASCE 109(1), Table 2)
c
       data cpwfdat /
     >      0.0495,     ! sand
     >      0.0613,     ! loamy sand
     >      0.1101,     ! sandy loam
     >      0.0889,     ! loam
     >      0.1668,     ! silt loam
     >      0.2185,     ! sandy clay loam
     >      0.2088,     ! clay loam
     >      0.2730,     ! silty clay loam
     >      0.2390,     ! sandy clay
     >      0.2922,     ! silty clay
     >      0.3163      ! clay
     >      /     
c
c -----------------------------------------------------------------------
c -----------------------------------------------------------------------
c added for Green-Ampt
c      do l = 1, npoi
c        fwpudtot(l) = 0.0
c        wpud(l) = 0.0
c      end do
c
c -----------------------------------------------------------------------
c set sand/silt/clay vectors (xdat,ydat,zdat) for 11 data points
c
      do 100 l = 1, ndat
        xdat(l) = texdat(1,l)
        ydat(l) = texdat(2,l)
        zdat(l) = texdat(3,l)
 100  continue
c
c initialization and normalization constant for puddle model (kg m-2)
c
c Changed by TET
c      wpudmax = 20.0
*      wpudmax = 10.0
c
      if (irestart .eq. 0) then
        call const (wpud,  npoi, 0.0)
        call const (wipud, npoi, 0.0)
      end if
c
c set prescribed soil layer thicknesses
c
*      hsoi(1)  = 0.05
*      hsoi(2)  = 0.05
*      hsoi(3)  = 0.1
*      hsoi(4)  = 0.1
*      hsoi(5)  = 0.1
*      hsoi(6)  = 0.2
*      hsoi(7)  = 0.2
*      hsoi(8)  = 0.2
*      hsoi(9)  = 0.5
*      hsoi(10) = 0.5
*      hsoi(11) = 0.5
c
c set physical parameters of soil
c
      call const (z0soi, npoi, 0.005)
c
c initialize soil water and soil temperature fields
c also initialize daily average fields for crops
c
      if (irestart .eq. 0) then
c
c
        call const (wsoi,  npoi*nsoilay, 0.50  )
        call const (wisoi, npoi*nsoilay, 0.00  )
        call const (tsoi,  npoi*nsoilay, 278.13)
c
        call const (adwsoilay,  npoi*nsoilay, 0.50  )
        call const (adwisoilay, npoi*nsoilay, 0.00  )
        call const (adtsoilay,  npoi*nsoilay, 278.13)
c
c initialize total irrigation
c
        call const (totirrig, npoi, 0.0)
c
c initialize plant available nitrogen pool (kg/m2)  
c
        call const (aplantn, npoi, 0.0080)
c
c initialize soil  nitrogen variables
c
        call const (totnuptake, npoi*npft, 0.0) 
        call const (stressn,    npoi*npft, 1.0) 
c
        call const (tg, npoi, 278.13)
        call const (ti, npoi, 273.13)
      else
         do 150 i = 1, npoi
            tg(i) = tsoi(i,1)
            ti(i) = tsno(i,1)
 150     continue
c
      end if
c
c set soil surface parameters for the global domain
c
c SEE BELOW OTHER PARAMETERS
c
c      do 200 i = 1, npoi
c
c Convert input sand and clay percents to fractions
c
c        msand = nint(sand(i,1))
c        mclay = nint(clay(i,1)) 
c
c        fsand = 0.01 * msand
c        fclay = 0.01 * mclay
c        fsilt = 0.01 * (100 - msand - mclay)
c
c soil surface albedo:
c
c from bats table 3.ii assuming albedo depends on texture
c
c        albsav(i) = fsand * 0.120 +
c     >              fsilt * 0.085 +
c     >              fclay * 0.050
c
c        albsan(i) = 2.0 * albsav(i)
c
c 200  continue   
c
c create soil properties look-up table
c
c set soil parameters at each layer for the global domain
c
!---------------------------note-------------------------------------------
!  when isoil_layer = 1, set domtexti = f(lon,lat,layer) from read-in NetCDF file
!  This is done by overwriting values in domtext(npoi,nsoilay)
!
!  Y.Li 2014-06-18
!
!-------------------------End note-----------------------------------------
!
      if(isoil_layer.eq.1) then
        !read-in NetCDF file
        call read_NetCDF('soil',0)
 
        !assign data to domtext
        do l =1, nsoilay
          call arr2vec(cdummy((l-1)*nlonsub*nlatsub + 1), domtext(:,l))
        end do

        if(flg_sandlayers == 1) then !modify last 2 layers to sand
          domtext(:,10) = 1
          domtext(:,11) = 1
        end if

      end if !on if(isoil_layer.eq.1)


      do 300 k = 1, nsoilay 
        do 310 i = 1, npoi 
c
          domtexti = nint(domtext(i,k))
c
c         domtexti = 12 
c
c set soil texture constant
c
c CONUS textural classes
c
c class number           class
c ------------           --------------
c   1                    sand
c   2                    loamy sand
c   3                    sandy loam
c   4                    silt loam
c   5                    silt
c   6                    loam
c   7                    sandy clay loam
c   8                    silty clay loam
c   9                    clay loam
c   10                   sandy clay
c   11                   silty clay
c   12                   clay
c 
c
c Convert input sand and clay percents to fractions
c
c Changed by TET
c          if (k.le.6) then
c          if (k.le.4) then
c            msand = nint(sand(i,k))
c            mclay = nint(clay(i,k)) 
c          else
c            msand = nint(sand(i,6)) 
c            mclay = nint(clay(i,6)) 
c            msand = nint(sand(i,4)) 
c            mclay = nint(clay(i,4)) 
c          endif
c
c          fsand    = 0.01 * msand
c          fclay    = 0.01 * mclay
c          fsilt    = 0.01 * (100 - msand - mclay)
c
c write out for each grid cell soil texture for layers
c
c          soisand(i,k) = fsand
c          soiclay(i,k) = fclay 
c
c for now, we assume that all soils have a 1% organic content -- 
c this is just a place holder until we couple the soil carbon
c dynamics to the soil physical properties
c
          forganic = 0.010
c
c density of soil material (without pores, not bulk) (kg m-3)
c from Campbell and Norman, 1998
c
          rhosoi(i,k) = 2650.0 * (1.0 - forganic) +
     >                  1300.0 * forganic 
c
c specific heat of soil material (j kg-1 k-1):
c from Campbell and Norman, 1998
c
          csoi(i,k) =  870.0 * (1.0 - forganic) +
     >                1920.0 * forganic 
c
c tet changed from cjk
c used to match textural fractions with soil textural class 
c calls two functions to match sand and clay fractions
c with proper soil textural class based on the usda
c classification system
c
c Now use input dominant soil texture class from soil_text.nc
c Must convert from CONUS texture class values to those of IBIS
c no silt in IBIS so set CONUS silt equal to silt loam
c
c          lmin = textcls (msand,mclay)
c

          if (domtexti .eq. 1) then
            lmin = 1
          else if (domtexti .eq. 2) then
            lmin = 2
          else if (domtexti .eq. 3) then
            lmin = 3
          else if (domtexti .eq. 4) then
            lmin = 5
          else if (domtexti .eq. 5) then
            lmin = 5
          else if (domtexti .eq. 6) then
            lmin = 4
          else if (domtexti .eq. 7) then
            lmin = 6
          else if (domtexti .eq. 8) then
            lmin = 8
          else if (domtexti .eq. 9) then
            lmin = 7
          else if (domtexti .eq. 10) then
            lmin = 9
          else if (domtexti .eq. 11) then
            lmin = 10
          else if (domtexti .eq. 12) then
            lmin = 11
          endif

          fracsand(i,k) = texdat(1,lmin)
          fracsilt(i,k) = texdat(2,lmin)
          fracclay(i,k) = texdat(3,lmin)
c porosity (fraction):
c 
          poros(i,k) = porosdat(lmin)
c
c field capacity (defined relative to the porosity):
c
          sfield(i,k) = 1.0 / poros(i,k) * sfielddat(lmin)
c
c wilting point (defined relative to the porosity):
c
          swilt(i,k)  = 1.0 / poros(i,k) * swiltdat(lmin)
c
c "b" exponent for the Campbell moisture-release equation:
c
          bex(i,k) = bexdat(lmin)
c
c nearest integer of "b" exponent (for computational efficiency):
c
          ibex(i,k) = nint(bex(i,k))
c
c saturated matric (air entry) potential (m-h2o):
c
          suction(i,k) = suctiondat(lmin)
c
c saturated hydraulic conductivity (m s-1):
c
          hydraul(i,k) = hydrauldat(lmin)
c
c
c -----------------------------------------------------------------------
c added for Green-Ampt
          cpwf(i, k) = cpwfdat(lmin)
          swater(i, k) = 0.000001   ! useless value
          sice(i, k) = 0        ! useless value
c -----------------------------------------------------------------------
c

 310    continue
 300  continue
c
c surface parameters
c
      do 200 i=1,npoi
        albsav(i) = fracsand(i,1) * 0.120 +
     >              fracsilt(i,1) * 0.085 +
     >              fracclay(i,1) * 0.050
c
        albsan(i) = 2.0 * albsav(i)
 200  continue
c
c return to main program
c
      return
      end
c
c
c ---------------------------------------------------------------------
      subroutine inimet (iyear, imonth, iday, nday, time, seed)
c ---------------------------------------------------------------------
c
c initializes surface meteorology from hourly data  
c avoids the need to go through much of subroutine 'daily' in weather.f 
c
      use comgrid
      use compar
      use comwork
      use comatm
      use comhyd
      use comsum
      use combcs
      use comhour
      use comcrop
      use comnitr
c
      implicit none
c
      integer nhours,
     >        npts,
     >        seed,
     >        nday,
     >        iyear,
     >        imonth,
     >        iday,
     >        ihtime,
     >        id1,
     >        id2,
     >        id3,  
     >        n,
     >        i,
     >        ihprecip,
     >        iyr,
     >        iyrinc
c      
      real time,
     >     rwork,
     >     var3sum,
     >     var3max,
     >     var3min,
     >     chtime,
     >     psurfi,
     >     tdave,
     >     shum,
     >     elevm, 
     >     ran2 
c
      character*100 header
      character*100 metfile
      character*30  hstg
      character*39  pathn   ! has to be exact length of input string
c
c External
c
      logical is_leap        ! Function: is the given year a leap year?
c
      save nhours,
     >     var3sum,
     >     var3max,
     >     var3min
c
      include 'comsat.h'
c
      pathn = '/Volumes/wet/kucharik/IBIS_input/input/'
c
c spin up period by choosing a random year
c of hourly weather from 1986-2000
c
c have to make sure that npts matches the number
c of hourly points in a particular year (iyear)
c during spin up
c
      npts = 8760  ! number of yourly points in non-leap year 
c
c generate the hourly met file to use (iyr)
c spin-up
c
      if (iyear .lt. 1986 .or. iyear .gt. 2000) then
c
c is the actual calendar year (iyear) a leap year?
c if so, change the number of points to read in
c
        if (is_leap(iyear)) then
          npts = 8784
        end if
c
c then choose a leap year out of met record  
c but only from 92, 96, or 2000 since 88 was
c an anomoly year (drought)
c
        if (npts .eq. 8784) then 
           iyrinc = min(8, int(ran2(seed) * 9))
           iyrinc = int(iyrinc/4) * 4
           iyr    = 1992 + iyrinc 
c
c not a calendar leap year - choose any year
c choose any other year - but note, if you were
c to choose a leap year, it will only read in the
c first 8760 pts from the met file
c
        else
          iyrinc = min(14, int(ran2(seed) * 15)) 
          iyrinc = max(0, iyrinc)
          iyr    = 1986 + iyrinc 
c
        endif
c
      else
c
c now met files match calendar year 
c
c change npts for leap years
c
        if (iyear .eq. 1988 .or.
     >      iyear .eq. 1992 .or.
     >      iyear .eq. 1996 .or.
     >      iyear .eq. 2000) npts = 8784             
c
        iyr = iyear
c
      endif
c

c arlington met station data files
c------------------------------------------
      if (iyr .eq. 1986) then 
         metfile = pathn//'arlington/metdata/1986'
      endif
c
      if (iyr .eq. 1987) then
         metfile = pathn//'arlington/metdata/1987'
      endif
c
      if (iyr .eq. 1988) then
         metfile = pathn//'arlington/metdata/1988'
      endif
c
      if (iyr .eq. 1989) then 
         metfile = pathn//'arlington/metdata/1989'
      endif
c
      if (iyr .eq. 1990) then 
         metfile = pathn//'arlington/metdata/1990'
      endif
c
      if (iyr .eq. 1991) then 
         metfile = pathn//'arlington/metdata/1991'
      endif
c
      if (iyr .eq. 1992) then
         metfile = pathn//'arlington/metdata/1992_lancaster'
      endif
c
      if (iyr .eq. 1993) then
         metfile = pathn//'arlington/metdata/1993'
      endif
c
      if (iyr .eq. 1994) then
         metfile = pathn//'arlington/metdata/1994'
      endif
c
      if (iyr .eq. 1995) then
         metfile = pathn//'arlington/metdata/1995'
      endif
c
      if (iyr .eq. 1996) then
         metfile = pathn//'arlington/metdata/1996'
      endif
c
      if (iyr .eq. 1997) then
         metfile = pathn//'arlington/metdata/1997'
      endif
c
      if (iyr .eq. 1998) then
         metfile = pathn//'arlington/metdata/1998'
      endif
c
      if (iyr .eq. 1999) then
         metfile = pathn//'arlington/metdata/1999'
      endif
c
      if (iyr .eq. 2000) then
         metfile = pathn//'arlington/metdata/2000'
      endif
c
c-----------------------------------------------
      if (imonth .eq. 1 .and. iday .eq. 1 .and.
     >    time .eq. 0.0) then
        write(6,*) 'Reading hourly weather data ...'
        write(6,*) 'Using met file year = ', iyr
        open(12,file = metfile)
        read(12,*) header
        do 100 n = 1, npts 
          read(12,500) iiyear(n),iimonth(n),iiday(n), iihour(n)
 100    continue
 500    format(i4,1x,i2,1x,i2,1x,i4)
c
        close(12)    
c
        open(12,file = metfile)
        read(12,*) header
        do 200 n = 1, npts 
          read(12,*) hstg,id1,var1(n),var2(n),var3(n),
     >                       var4(n),var5(n),var6(n)
 200    continue
        close(12)
c
      endif
c
        rwork = (grav / rair / 0.0065)
c
c set elevation for the arlington grid cells
c based on elevation reading (meters) at morrisonville
c
        elevm = 293.8
c
        do 400 n  = 1, npts 
c
c         if (iyear .eq. iiyear(n)) then
            if (imonth .eq. iimonth(n)) then
              if (iday .eq. iiday(n))     then
c
                 chtime = (time / 3600.)
                 ihtime = int(chtime) + 1
c
                if (ihtime .eq. (iihour(n) / 100)) then
                  if (time .eq. 0.) then
                    var3sum = 0.0
                    nhours  = 0
                    call const (tmin, npoi, 400.0)
                    call const (tmax, npoi, 200.0)
                  endif              
c
                  nhours  = nhours + 1
                  var3min = min(tmin(1), var3(n) + 273.16)
                  var3max = max(tmax(1), var3(n) + 273.16)
                  psurfi  = 101325.0 * ((var3(n) + 273.16) /
     >                      ((var3(n) + 273.16) + 0.0065 * 
     >                        elevm)) ** rwork

                  call const (psurf, npoi, psurfi)         ! surface pressure
c
c for precipitation, value at time n is precip accumulated since last
c time step, so we need to get the precip value for the time step n+1 
c
                  ihprecip = min(n+1, npts)  ! use precip from the next time step
                  call const (precip, npoi, var1(ihprecip)
     >                        * dtime / 3600.)             ! precipitation (mm)  
                  call const (tmin, npoi, var3min)         ! minimum temp
                  call const (tmax, npoi, var3max)         ! maximum temp
                  call const (ta, npoi, var3(n) + 273.16)  ! air temperature
                  var3sum = var3sum + var3(n) + 273.16
                  call const (cloud, npoi, var2(n))        ! insolation
c
                  if (var4(n) .lt. 1.) then
                    print *, 'rh too low! ', var4(n)
                    print *, 'month, day, hour ', imonth, iday, ihtime,n
                    stop 33
                  endif 
c               
c convert relative humidity to specific humidity
c
                  
                  shum  = qsat((var4(n)/100.)*esat(var3(n) + 273.16),
     >                      psurfi) 
c
                  call const (qd, npoi, shum)             ! specific humidity 
                  call const (ua, npoi, var5(n))           ! wind velocity
c
                  if (time .eq. 24. * 3600. - dtime) then
                    if (nhours .ne. 24. * int (3600./dtime)) then
                      print *, 'this day has incorrect hours!'
                      print *, 'nhours = ', nhours
                      print *, 'iyear ', iyr 
                      stop 7
                    endif
                    var3sum = var3sum / float(nhours)
c
c three ways to compute daily average temp
c * average of hours        :  tdave = var3sum
c * weighted of min and max :  tdave = tmax_wt * var3max + tmin_wt * var3min
c * average of min and max  :  tdave = (var3min + var3max) / 2 
c
                   tdave = var3sum
                   call const (td, npoi, tdave)
c                  call dailymet(imonth, iday, seed, tdave, 0)
                   call dailymet(imonth, iday, seed, 1)
c
                 endif     
                endif
              endif
            endif
c         endif
c
 400    continue
c
c
c return to main program
c
      return
      end
c
c-------------------------------------------------------------------------
      integer function textcls (msand,mclay)
c
c adapted for ibis by cjk 01/11/01
c-------------------------------------------------------------------------
c |
c |                         T R I A N G L E
c | Main program that calls WHAT_TEXTURE, a function that classifies soil
c | in the USDA textural triangle using sand and clay %
c +-----------------------------------------------------------------------
c | Created by: aris gerakis, apr. 98 with help from brian baer
c | Modified by: aris gerakis, july 99: now all borderline cases are valid
c | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
c |              main program
c +-----------------------------------------------------------------------
c | COMMENTS
c | o Supply a data file with two columns, in free format:  1st column sand,
c |   2nd column clay %, no header.  The output is a file with the classes.
c +-----------------------------------------------------------------------
c | You may use, distribute and modify this code provided you maintain
c ! this header and give appropriate credit.
c +-----------------------------------------------------------------------
c
c code adapted for IBIS by cjk 01-11-01
c
c     use comgrid
c     use compar
c     use comsoi
c
      integer  msand,
     >         mclay
c
      logical inpoly
c
      real    silty_loam(1:7,1:2),
     >        sandy(1:7,1:2),
     >        silty_clay_loam(1:7,1:2), 
     >        loam(1:7,1:2),
     >        clay_loam(1:7,1:2),
     >        sandy_loam(1:7,1:2),
     >        silty_clay(1:7,1:2),
     >        sandy_clay_loam(1:7,1:2), 
     >        loamy_sand(1:7,1:2),
     >        clayey(1:7,1:2),
c    >        silt(1:7,1:2), 
     >        sandy_clay(1:7,1:2)
c
c initalize polygon coordinates:
c each textural class reads in the sand coordinates (1,7) first, and
c then the corresponding clay coordinates (1,7)

c     data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
c
c because we do not have a separate silt category, have to redefine the
c polygon boundaries for the silt loam  
c
      data sandy           /85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
      data loamy_sand      /70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
      data sandy_loam      /50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 15, 0/
      data loam            /43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
      data silty_loam      /0, 0, 23, 50, 0, 0, 0, 0, 27, 27, 0, 0, 0, 0/ 
c     data silt            /0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
      data sandy_clay_loam /52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35, 20, 0, 0/
      data clay_loam       /20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data silty_clay_loam /0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data sandy_clay      /45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/
      data silty_clay      /0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
      data clayey          /20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/
c
c polygon coordinates  
c
c     sand
c
c     >  85, 90, 100, 0, 0, 0, 0,       ! sand
c     >  70, 85, 90, 85, 0, 0, 0,       ! loamy sand
c     >  50, 43, 52, 52, 80, 85, 70,    ! sandy loam
c     >  43, 23, 45, 52, 52, 0, 0,      ! loam
c     >   0, 0, 23, 50, 0, 0, 0,        ! silt loam (combined with silt)
c     >  52, 45, 45, 65, 80, 0, 0,      ! sandy clay loam
c     >  20, 20, 45, 45, 0, 0, 0,       ! clay loam
c     >   0, 0, 20, 20, 0, 0, 0,        ! silty clay loam
c     >  45, 45, 65, 0, 0, 0, 0,        ! sandy clay
c     >   0, 0, 20, 0, 0, 0, 0,         ! silty clay 
c     >  20, 0, 0, 45, 45, 0, 0         ! clay
c
c      clay
c
c     > 0, 10, 0, 0, 0, 0, 0,           ! sand
c     > 0, 15, 10, 0, 0, 0, 0,          ! loamy sand
c     > 0, 7, 7, 20, 20, 15, 0,         ! sandy loam 
c     > 7, 27, 27, 20, 7, 0, 0,         ! loam
c     > 0, 27, 27, 0, 0, 0, 0,          ! silt loam (combined with silt)
c     > 20, 27, 35, 35, 20, 0, 0,       ! sandy clay loam
c     > 27, 40, 40, 27, 0, 0, 0,        ! clay loam
c     > 27, 40, 40, 27, 0, 0, 0,        ! silty clay loam
c     > 35, 55, 35, 0, 0, 0, 0,         ! sandy clay
c     > 40, 60, 40, 0, 0, 0, 0,         ! silty clay
c     > 40, 60, 100, 55, 40, 0, 0       ! clay
c
c +-----------------------------------------------------------------------
c | figure out what texture grid cell and layer are part of  
c | classify a soil in the triangle based on sand and clay %
c +-----------------------------------------------------------------------
c | Created by: aris gerakis, apr. 98
c | Modified by: aris gerakis, june 99.  Now check all polygons instead of
c | stopping when a right solution is found.  This to cover all borderline 
c | cases.
c +-----------------------------------------------------------------------
c
c C.Molling 8/18/2004
c The following lines of code are necessary, as the algorithm always
c gives siltloam if sand or clay is 0%.  We use a minimum of 1% for
c sand and clay, and then adjust the greater of the two if the sum
c is more than 100%.
c
      msand = max(msand,1)
      mclay = max(mclay,1)
      if (msand + mclay .gt. 100) then
         if (msand .gt. mclay) then
            msand = msand - 1
         else
            mclay = mclay - 1
         endif
      endif
c
c end C. Molling 8/18/2004
c
c find polygon(s) where the point is.  
c
      textcls = 0 
c
      if (msand .gt. 0.0 .and. mclay .gt. 0.0) then
         if (inpoly(sandy, 3, msand, mclay)) then
            textcls = 1      ! sand
         endif
         if (inpoly(loamy_sand, 4, msand, mclay)) then
            textcls = 2      ! loamy sand
         endif
         if (inpoly(sandy_loam, 7, msand, mclay)) then
            textcls = 3      ! sandy loam
         endif
         if (inpoly(loam, 5, msand, mclay)) then
            textcls = 4      ! loam
         endif
         if (inpoly(silty_loam, 4, msand, mclay)) then
            textcls = 5      ! silt loam
         endif
         if (inpoly(sandy_clay_loam, 5, msand, mclay)) then
            textcls = 6      ! sandy clay loam
         endif
         if (inpoly(clay_loam, 4, msand, mclay)) then
            textcls = 7      ! clay loam
         endif
         if (inpoly(silty_clay_loam, 4, msand, mclay)) then
            textcls = 8      ! silty clay loam
         endif
         if (inpoly(sandy_clay, 3, msand, mclay)) then
            textcls = 9      ! sandy clay
         endif
         if (inpoly(silty_clay, 3, msand, mclay)) then
            textcls = 10     ! silty clay
         endif
         if (inpoly(clayey, 5, msand, mclay)) then
            textcls = 11     ! clay
         endif
      endif
c
      if (textcls .eq. 0) then
         textcls = 5         ! silt loam
c
c        write (*, 1000) msand, mclay
c 1000   format (/, 1x, 'Texture not found for ', f5.1, ' sand and ', f5.1, ' clay')
      endif
c
      return
      end
c
c---------------------------------------------------------------------------
      logical function inpoly (poly, npoints, xt, yt)
c
c adapted for ibis by cjk 01/11/01
c---------------------------------------------------------------------------
c
c                             INPOLY
c !   Function to tell if a point is inside a polygon or not.
c !--------------------------------------------------------------------------
c !   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
c !
c !   Please feel free to use this source code for any purpose, commercial
c !   or otherwise, as long as you don't restrict anyone else's use of
c !   this source code.  Please give credit where credit is due.
c !
c !   Point-in-polygon algorithm, created especially for World-Wide Web
c !   servers to process image maps with mouse-clickable regions.
c !
c !   Home for this file:  http://www.gcomm.com/develop/inpoly.c
c !
c !                                       6/19/95 - Bob Stein & Craig Yap
c !                                       stein@gcomm.com
c !                                       craig@cse.fau.edu
c !--------------------------------------------------------------------------
c !   Modified by:
c !   Aris Gerakis, apr. 1998: 1.  translated to Fortran
c !                            2.  made it work with real coordinates
c !                            3.  now resolves the case where point falls
c !                                on polygon border.
c !   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
c !   Aris Gerakis, july 1999: Now all borderline cases are valid
c !--------------------------------------------------------------------------
c !   Glossary:
c !   function inpoly: true=inside, false=outside (is target point inside
c !                    a 2D polygon?)
c !   poly(*,2):  polygon points, [0]=x, [1]=y
c !   npoints: number of points in polygon
c !   xt: x (horizontal) of target point
c !   yt: y (vertical) of target point
c !--------------------------------------------------------------------------
c
c declare arguments  
c
      integer  npoints,
     >         xt,
     >         yt 
c
      real     poly(7, 2)
c
c local variables
c
      real     xnew,
     >         ynew,
     >         xold,
     >         yold,
     >         x1,
     >         y1,
     >         x2,
     >         y2
c
      integer  i
c
      logical inside,
     >        on_border

      inside = .false.
      on_border = .false.
c
      if (npoints .lt. 3)  then
        inpoly = .false.
        return
      end if
c
      xold = poly(npoints,1)
      yold = poly(npoints,2)

      do 300  i = 1 , npoints
        xnew = poly(i,1)
        ynew = poly(i,2)

        if (xnew .gt. xold)  then
          x1 = xold
          x2 = xnew
          y1 = yold
          y2 = ynew
        else
          x1 = xnew
          x2 = xold
          y1 = ynew
          y2 = yold
        end if

c the outer IF is the 'straddle' test and the 'vertical border' test.
c the inner IF is the 'non-vertical border' test and the 'north' test.  

c the first statement checks whether a north pointing vector crosses  
c (stradles) the straight segment.  There are two possibilities, depe-
c nding on whether xnew < xold or xnew > xold.  The '<' is because edge 
c must be "open" at left, which is necessary to keep correct count when 
c vector 'licks' a vertix of a polygon.  

        if ((xnew .lt. xt .and. xt .le. xold) .or. (.not. xnew .lt. xt .and. 
     >     .not. xt .le. xold)) then
c
c the test point lies on a non-vertical border:
c
          if ((yt-y1)*(x2-x1) .eq. (y2-y1)*(xt-x1)) then
            on_border = .true. 
c
c check if segment is north of test point.  If yes, reverse the 
c value of INSIDE.  The +0.001 was necessary to avoid errors due   
c arithmetic (e.g., when clay = 98.87 and sand = 1.13):   
c
          elseif ((yt-y1)*(x2-x1) .lt. (y2-y1)*(xt-x1) + 0.001) then
            inside = .not.inside ! cross a segment
          endif
c
c this is the rare case when test point falls on vertical border or  
c left edge of non-vertical border. The left x-coordinate must be  
c common.  The slope requirement must be met, but also point must be
c between the lower and upper y-coordinate of border segment.  There 
c are two possibilities,  depending on whether ynew < yold or ynew > 
c yold:
c
        elseif ((xnew .eq. xt .or. xold .eq. xt) .and. (yt-y1)*(x2-x1) .eq. 
     >    (y2-y1)*(xt-x1) .and. ((ynew .le. yt .and. yt .le. yold) .or. 
     >    (.not. ynew .lt. yt .and. .not. yt .lt. yold))) then
          on_border = .true. 
        endif
c
        xold = xnew
        yold = ynew
c
 300    continue  
c
c If test point is not on a border, the function result is the last state 
c of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
c inside the polygon if it falls on any of its borders:
c
      if (.not. on_border) then
         inpoly = inside
      else
         inpoly = .true.
      endif
c
      return
      end
c
c
c
c ---------------------------------------------------------------------
      subroutine iniveg (isimveg, irestart)
c ---------------------------------------------------------------------
c
      use comgrid
      use compar
      use combcs
      use comsoi
      use comsum
      use comveg
      use comcrop
      use comnitr
      use compft
c
      implicit none
c
c Arguments (input)
c
      integer irestart,     ! 0: not a restart run 1: restart run
     >        isimveg       ! 0:static veg 1:dyn veg 2:dyn veg-cold start
c
c local variables
c
      integer ideci,        ! # deciduous plant functional types (pft)
     >        ievgr,        ! # evergreen pft 
     >        ishrub,       ! # shrub pft 
     >        igrass,       ! # herbaceous pft
     >        icrop,        ! # cropland pft 
     >        ilower,       ! possible # pft for lower canopy
     >        iupper,       ! possible # pft for upper canopy
     >        inveg,        ! vegetation type
     >        i,j,k         ! loop indices
c
      real
     >     plaievgr,        ! potential lai of evergreen trees
     >     plaideci,        ! potential lai of deciduous trees
     >     plaishrub,       ! potential lai of shrubs
     >     plaigrass,       ! potential lai of grasses
*     >     plaiupper,       ! potential lai of upper canopy (uniform initial veg)
*     >     plailower,       ! potential lai of lower canopy (uniform initial veg)
*     >     xminlai,         ! minimum lai for each existing plant type
     >     wood,            ! total wood biomas in grid cell
     >     chiflz,          ! lower canopy leaf orientation factor 
     >                      !    (-1 vertical, 0 random, 1 horizontal)
     >     chifuz,          ! uppuer canopy leaf orientation factor 
*    >     beta1,           ! parameter for Jackson rooting profile, lower canopy
*    >     beta2,           ! parameter for Jackson rooting profile, upper canopy
     >     totdepth,        ! total soil depth
     >     frootnorm1,      ! normalization factor for Jackson rooting profile,low
     >     frootnorm2       ! normalization factor for Jackson rooting profile, up
     
c
      real depth(nsoilay)   ! soil layer depth (cm)
c
c initialize specific leaf area values
c
*      data specla  / 25.0,  ! tropical broadleaf evergreen trees
*     >               25.0,  ! tropical broadleaf drought-deciduous trees
*     >               25.0,  ! warm-temperate broadleaf evergreen trees
*     >               12.5,  ! temperate conifer evergreen trees
*     >               25.0,  ! temperate broadleaf cold-deciduous trees
*     >               12.5,  ! boreal conifer evergreen trees
*     >               25.0,  ! boreal broadleaf cold-deciduous trees  
*     >               25.0,  ! boreal conifer cold-deciduous trees
*     >               12.5,  ! evergreen shrubs 
*     >               25.0,  ! deciduous shrubs 
*     >               20.0,  ! warm (c4) grasses
*     >               20.0,  ! cool (c3) grasses
*     >               55.0,  ! soybean
*     >               25.0,  ! maize 
*     >               55.0 / ! wheat
c
c initial base temperature (K) values for gdd accumulation for crops
c
*      baset(13) =  283.16  ! soybean 
*      baset(14) =  281.16  ! maize 
*      baset(15) =  273.16  ! wheat 
c
c initial maximum daily gdd (C) accumulation for crops
c
*      mxtmp(13) =  30.0  ! soybean 
*      mxtmp(14) =  30.0  ! maize 
*      mxtmp(15) =  26.0  ! wheat 
c
*      woodnorm = 7.5
c
      if (irestart .eq. 0) then
        call const     (cntops,     npoi*npft, 40.0) 
        call const     (cnroot,     npoi*npft, 60.0)
        call const     (tnplant,    npoi*npft, 0.0) 
        call const_int (harvdate,   npoi*npft, 999)
      endif
c
      do 100 i = 1, npoi
        if(cmask(i)==0) cycle
c
c initialize the moisture stress factors
c
        stresstu(i) = 1.0
        stresstl(i) = 1.0
c
c initialize running-mean air temperature
c
        if (irestart .eq. 0) then
c
           a10td(i) = 273.16
           a11soiltd(i) = 273.16
           daylength(i) = 0.
           a3tdmin(i) = 0.
           precipsum(i) = 0.
           stsumu(i) = 0.
           stsuml(i) = 0.

c
c initialize running-mean values of canopy photosynthesis rates
c
           a10ancub(i) = 10.0e-06
           a10ancuc(i) = 10.0e-06
           a10ancls(i) = 10.0e-06
           a10ancl4(i) = 10.0e-06
           a10ancl3(i) = 10.0e-06
c 
c initialize running-mean values of the scaling parameter
c
           a10scalparamu(i) = 0.5 * 5.
           a10scalparaml(i) = 0.5 * 5.
           a10daylightu(i) = 5.
           a10daylightl(i) = 5.
c
c initialize litter fall
c
           falll(i) = 0.
           fallr(i) = 0.
           fallw(i) = 0.
c
        end if
c
c reset counters
c
        ievgr  = 0
        ideci  = 0
        ishrub = 0
        igrass = 0
        icrop  = 0
        ilower = 0
        iupper = 0
c
c determine number of evergreen plant functional types
c
        if (nint(exist(i,1)).eq.1) ievgr = ievgr + 1
        if (nint(exist(i,3)).eq.1) ievgr = ievgr + 1
        if (nint(exist(i,4)).eq.1) ievgr = ievgr + 1
        if (nint(exist(i,6)).eq.1) ievgr = ievgr + 1
c
c determine number of deciduous plant functional types
c
        if (nint(exist(i,2)).eq.1) ideci = ideci + 1
        if (nint(exist(i,5)).eq.1) ideci = ideci + 1
        if (nint(exist(i,7)).eq.1) ideci = ideci + 1
        if (nint(exist(i,8)).eq.1) ideci = ideci + 1
c
c make sure counter is at least 1 (to avoid division by zero)
c
        ievgr = max (1, ievgr)
        ideci = max (1, ideci)
c
c determine number of shrub functional types
c
        if (nint(exist(i,9)).eq.1)  ishrub = ishrub + 1
        if (nint(exist(i,10)).eq.1) ishrub = ishrub + 1
c
c determine number of herbaceous plant functional types
c
        if (nint(exist(i,11)).eq.1) igrass = igrass + 1
        if (nint(exist(i,12)).eq.1) igrass = igrass + 1
c
c make sure counter is at least 1 (to avoid division by zero)
c
        ishrub = max (1, ishrub)
        igrass = max (1, igrass)
c
c total number of possible pfts for each canopy
c
        iupper = ievgr  + ideci
        ilower = ishrub + igrass 
c
c make sure counter is at least 1 (to avoid division by zero)
c
c determine number of crop functional types
c
        if (nint(exist(i,13)).eq.1) then
            icrop = icrop + 1   ! soybean 
            croplive(i,13) = 0.0 
            croppresent(i,13) = 0
        endif
        if (nint(exist(i,14)).eq.1) then  
            icrop = icrop + 1   ! maize 
            croplive(i,14) = 0.0 
            croppresent(i,14) = 0
        endif
        if (nint(exist(i,15)).eq.1) then 
            icrop = icrop + 1   ! wheat 
            croplive(i,15) = 0.0 
            croppresent(i,15) = 0
        endif
c
c make sure counter is at least 1 (to avoid division by zero)
c
        icrop = max (1, icrop)
        ilower = ilower + icrop
c
        iupper = max (1, iupper)
        ilower = max (1, ilower)
c
c cold start of the vegetation
c
        if (irestart .eq. 0) then
c
c
c ************************************************************************
c case (0) assign vegetation characteristics for static vegetation
c ************************************************************************
c
c and
c
c ************************************************************************
c case (1) assign vegetation characteristics for dynamic vegtation
c          that is initialized with fixed vegetation map
c ************************************************************************
c
          if (isimveg.eq.0 .or. isimveg.eq.1) then
c
c translate vegetation type (real) to nearest integer
c
             inveg = nint (xinveg(i))
c
c for initialization purposes, set the predicted vegetation type
c to the initial vegetation type
c
            vegtype0(i) = xinveg(i)
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
c these classes consist of some combination of 
c plant functional types:
c
c ---------------------------------------------------
c  1: tropical broadleaf evergreen trees
c  2: tropical broadleaf drought-deciduous trees
c  3: warm-temperate broadleaf evergreen trees
c  4: temperate conifer evergreen trees
c  5: temperate broadleaf cold-deciduous trees
c  6: boreal conifer evergreen trees
c  7: boreal broadleaf cold-deciduous trees
c  8: boreal conifer cold-deciduous trees
c  9: evergreen shrubs
c 10: cold-deciduous shrubs
c 11: warm (c4) grasses
c 12: cool (c3) grasses
c 13: soybean
c 14: maize 
c 15: winter and spring wheat 
c ---------------------------------------------------
c
c initially all values are set to zero
c

**** DTP 2001/05/25. The following code replaces the 450+
*    lines of stuff that follows it (hence the temporary goto
*    statement). Note that values of plai_init are read in as
*    parameters from params.veg. Note also that the declarations
*    of the four local variables plaievgr, plaideci, plaishrub 
*    and plaigrass can all be dropped.

            if(inveg/=0) then !Y.Li

              plai(i,1)  = exist(i,1)  / float(ievgr)  * 
     >                     plai_init(1,inveg)
              plai(i,2)  = exist(i,2)  / float(ideci)  * 
     >                     plai_init(2,inveg)
              plai(i,3)  = exist(i,3)  / float(ievgr)  * 
     >                     plai_init(1,inveg)
              plai(i,4)  = exist(i,4)  / float(ievgr)  * 
     >                     plai_init(1,inveg)
              plai(i,5)  = exist(i,5)  / float(ideci)  * 
     >                     plai_init(2,inveg)
              plai(i,6)  = exist(i,6)  / float(ievgr)  * 
     >                     plai_init(1,inveg)
              plai(i,7)  = exist(i,7)  / float(ideci)  * 
     >                     plai_init(2,inveg)
              plai(i,8)  = exist(i,8)  / float(ideci)  * 
     >                     plai_init(2,inveg)
              plai(i,9)  = exist(i,9)  / float(ishrub) * 
     >                     plai_init(3,inveg)
              plai(i,10) = exist(i,10) / float(ishrub) * 
     >                     plai_init(3,inveg)
*              plai(i,11) = exist(i,11) / float(igrass) * 
*    >                     plai_init(4,inveg)
*              plai(i,12) = exist(i,12) / float(igrass) * 
*     >                     plai_init(4,inveg)
              if ((inveg.eq.9).or.(inveg.eq.10)) then
                if (tw(i).gt.22.0) then
                  plai(i,11) = exist(i,11) * 0.80 * plai_init(4,inveg)
                  plai(i,12) = exist(i,12) * 0.20 * plai_init(4,inveg)
                else
                  plai(i,11) = exist(i,11) * 0.00 * plai_init(4,inveg)
                  plai(i,12) = exist(i,12) * 1.00 * plai_init(4,inveg)
                endif
              else
                plai(i,11) = exist(i,11) / float(igrass) *
     >                     plai_init(4,inveg)
                plai(i,12) = exist(i,12) / float(igrass) *
     >                     plai_init(4,inveg)
              endif

            end if  !on if(inveg/=0)

* cropping systems - just set equal to zero
*
            plai(i,13) = 0.0  ! soybean
            plai(i,14) = 0.0  ! maize
            plai(i,15) = 0.0  ! wheat

      GOTO 9999


            plai(i,1)  = 0.0
            plai(i,2)  = 0.0
            plai(i,3)  = 0.0
            plai(i,4)  = 0.0
            plai(i,5)  = 0.0
            plai(i,6)  = 0.0
            plai(i,7)  = 0.0
            plai(i,8)  = 0.0
            plai(i,9)  = 0.0
            plai(i,10) = 0.0
            plai(i,11) = 0.0
            plai(i,12) = 0.0
            plai(i,13) = 0.0
            plai(i,14) = 0.0
            plai(i,15) = 0.0
c
c ---------------------------------------------------
c  1: tropical evergreen forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.1) then
c
              plaievgr  = 5.00
              plaideci  = 1.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  2: tropical deciduous forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.2) then
c
              plaievgr  = 1.00
              plaideci  = 5.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  3: temperate evergreen broadleaf forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.3) then
c
              plaievgr  = 4.00
              plaideci  = 1.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  4: temperate evergreen conifer forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.4) then
c
              plaievgr  = 3.00
              plaideci  = 1.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  5: temperate deciduous forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.5) then
c
              plaievgr  = 1.00
              plaideci  = 3.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  6: boreal evergreen forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.6) then
c
              plaievgr  = 3.00
              plaideci  = 1.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  7: boreal deciduous forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.7) then
c
              plaievgr  = 1.00
              plaideci  = 3.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  8: mixed forest / woodland
c ---------------------------------------------------
c
            if (inveg.eq.8) then
c
              plaievgr  = 2.00
              plaideci  = 2.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c  9: savanna
c ---------------------------------------------------
c
            if (inveg.eq.9) then
c
              plaievgr  = 0.50
              plaideci  = 1.00
              plaishrub = 0.50
              plaigrass = 2.00
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
c
c determine if c4/c3 dominated grassland
c
              if (tw(i).gt.22.0) then
                plai(i,11) = exist(i,11) * 0.80 * plaigrass
                plai(i,12) = exist(i,12) * 0.20 * plaigrass
              else
                plai(i,11) = exist(i,11) * 0.00 * plaigrass
                plai(i,12) = exist(i,12) * 1.00 * plaigrass
              endif
c
            endif
c
c ---------------------------------------------------
c 10: grassland / steppe
c ---------------------------------------------------
c
            if (inveg.eq.10) then
c
              plaievgr  = 0.25
              plaideci  = 0.25
              plaishrub = 0.50
              plaigrass = 2.50
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
c
c determine if c4/c3 dominated savanna
c
              if (tw(i).gt.22.0) then
                plai(i,11) = exist(i,11) * 0.80 * plaigrass
                plai(i,12) = exist(i,12) * 0.20 * plaigrass
              else
                plai(i,11) = exist(i,11) * 0.00 * plaigrass
                plai(i,12) = exist(i,12) * 1.00 * plaigrass
              endif
c
            endif
c
c ---------------------------------------------------
c 11: dense shrubland
c ---------------------------------------------------
c
            if (inveg.eq.11) then
c
              plaievgr  = 0.10
              plaideci  = 0.10
              plaishrub = 1.00
              plaigrass = 0.50
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c 12: open shrubland
c ---------------------------------------------------
c
            if (inveg.eq.12) then
c
              plaievgr  = 0.00
              plaideci  = 0.00
              plaishrub = 0.25
              plaigrass = 0.25
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c 13: tundra
c ---------------------------------------------------
c
            if (inveg.eq.13) then
c
              plaievgr  = 0.00
              plaideci  = 0.00
              plaishrub = 1.00
              plaigrass = 1.00
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c 14: desert
c ---------------------------------------------------
c
            if (inveg.eq.14) then
c
              plaievgr  = 0.00
              plaideci  = 0.00
              plaishrub = 0.05
              plaigrass = 0.05
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
c ---------------------------------------------------
c 15: polar desert / rock / ice
c ---------------------------------------------------
c
            if (inveg.eq.15) then
c
              plaievgr  = 0.00
              plaideci  = 0.00
              plaishrub = 0.05
              plaigrass = 0.05
c
              plai(i,1)  = exist(i,1)  / float(ievgr)  * plaievgr
              plai(i,2)  = exist(i,2)  / float(ideci)  * plaideci
              plai(i,3)  = exist(i,3)  / float(ievgr)  * plaievgr
              plai(i,4)  = exist(i,4)  / float(ievgr)  * plaievgr
              plai(i,5)  = exist(i,5)  / float(ideci)  * plaideci
              plai(i,6)  = exist(i,6)  / float(ievgr)  * plaievgr
              plai(i,7)  = exist(i,7)  / float(ideci)  * plaideci
              plai(i,8)  = exist(i,8)  / float(ideci)  * plaideci
              plai(i,9)  = exist(i,9)  / float(ishrub) * plaishrub
              plai(i,10) = exist(i,10) / float(ishrub) * plaishrub
              plai(i,11) = exist(i,11) / float(igrass) * plaigrass
              plai(i,12) = exist(i,12) / float(igrass) * plaigrass
c
            endif
c
 9999      endif
c
c ************************************************************************
c case (2) assign vegetation characteristics for dynamic vegtation
c          that is initialized with uniform vegetation conditions
c ************************************************************************
c
c specify uniform initial conditions
c
          if (isimveg.eq.2) then
c
*            plaiupper = 0.5
*            plailower = 0.5
c
            plai(i,1)  = exist(i,1)  / float(iupper) * plaiupper
            plai(i,2)  = exist(i,2)  / float(iupper) * plaiupper
            plai(i,3)  = exist(i,3)  / float(iupper) * plaiupper
            plai(i,4)  = exist(i,4)  / float(iupper) * plaiupper
            plai(i,5)  = exist(i,5)  / float(iupper) * plaiupper
            plai(i,6)  = exist(i,6)  / float(iupper) * plaiupper
            plai(i,7)  = exist(i,7)  / float(iupper) * plaiupper
            plai(i,8)  = exist(i,8)  / float(iupper) * plaiupper
            plai(i,9)  = exist(i,9)  / float(ilower) * plailower
            plai(i,10) = exist(i,10) / float(ilower) * plailower
            plai(i,11) = exist(i,11) / float(ilower) * plailower
            plai(i,12) = exist(i,12) / float(ilower) * plailower
            plai(i,13) = exist(i,13) / float(ilower) * plailower
            plai(i,14) = exist(i,14) / float(ilower) * plailower
            plai(i,15) = exist(i,15) / float(ilower) * plailower

c
          endif
c
c ************************************************************************
c for both cases (1) and (2)
c ************************************************************************
c
c set minimum lai for each existing plant type
c
*         xminlai = 0.010
c
          plai(i,1)  = max (plai(i,1) , exist(i,1)  * xminlai)
          plai(i,2)  = max (plai(i,2) , exist(i,2)  * xminlai)
          plai(i,3)  = max (plai(i,3) , exist(i,3)  * xminlai)
          plai(i,4)  = max (plai(i,4) , exist(i,4)  * xminlai)
          plai(i,5)  = max (plai(i,5) , exist(i,5)  * xminlai)
          plai(i,6)  = max (plai(i,6) , exist(i,6)  * xminlai)
          plai(i,7)  = max (plai(i,7) , exist(i,7)  * xminlai)
          plai(i,8)  = max (plai(i,8) , exist(i,8)  * xminlai)
          plai(i,9)  = max (plai(i,9) , exist(i,9)  * xminlai)
          plai(i,10) = max (plai(i,10), exist(i,10) * xminlai)
          plai(i,11) = max (plai(i,11), exist(i,11) * xminlai)
          plai(i,12) = max (plai(i,12), exist(i,12) * xminlai)
          plai(i,13) = max (plai(i,13), exist(i,13) * xminlai)
          plai(i,14) = max (plai(i,14), exist(i,14) * xminlai)
          plai(i,15) = max (plai(i,15), exist(i,15) * xminlai)
c
c
c set sapwood fraction and biomass characteristics
c
*         sapfrac(i) = 0.1
          sapfrac(i) = sapfrac_init   ! from params.veg
c
          wood = 0.0
c
          do 120 j = 1, npft
c
c Note that, for croplands, these variables and other allocation
c fractions are re-initialized in initialcrop.
c
            aleaf(i,j) = aleaf_init(j)
            aroot(i,j) = aroot_init(j)
            awood(i,j) = awood_init(j)
c            
            cbiol(i,j) = plai(i,j) / specla(j)
            cbior(i,j) = 0.5 * cbiol(i,j)
c
            cbiow(i,j) = 0.0
c
c crop biomass storage -- stem and grain (fruit)
c
            cbios(i,j) = 0.0
            cbiog(i,j) = 0.0
c
c
            if (j.lt.9) cbiow(i,j) = plai(i,j) * 10.0 / 6.0
c
            biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j) +
     >                     cbios(i,j) + cbiog(i,j)
            wood = wood + cbiow(i,j)
c
 120      continue
c
c ************************************************************************
c case (3) assign vegetation characteristics for model restart
c ************************************************************************
c
        else    ! else for restart if loop
c
          wood = 0.0
c
          do 130 j = 1, npft
c
c WJS 06.03.10: Ideally, these allocation fractions would be saved to
c and read from the restart file, to allow for variations in time (as in
c croplands); however, as of now, they are not, so we need to initialize
c them to their initial values here.
c
c Note that, for croplands, these variables and other allocation
c fractions are re-initialized in initialcrop.
c
            aleaf(i,j) = aleaf_init(j)
            aroot(i,j) = aroot_init(j)
            awood(i,j) = awood_init(j)
c
            plai(i,j) = cbiol(i,j) * specla(j)
            biomass(i,j) = cbiol(i,j) + cbiow(i,j) + cbior(i,j)
            wood = wood + cbiow(i,j)
c
 130      continue
c
        end if  ! end restart if loop
c
c ************************************************************************
c determine basic vegetation structure characteristics
c ************************************************************************
c 
c total leaf area for upper and lower canopies
c
        totlaiu(i)  =  plai(i,1) + plai(i,2) +
     >                 plai(i,3) + plai(i,4) +
     >                 plai(i,5) + plai(i,6) +
     >                 plai(i,7) + plai(i,8) 
c
        totlail(i)  =  plai(i,9)  + plai(i,10) +
     >                 plai(i,11) + plai(i,12) +
     >                 plai(i,13) + plai(i,14) 
c
        totbiou(i)  = biomass(i,1) +
     >                biomass(i,2) +
     >                biomass(i,3) +
     >                biomass(i,4) +
     >                biomass(i,5) +
     >                biomass(i,6) +
     >                biomass(i,7) +
     >                biomass(i,8) 
c
        totbiol(i)  = biomass(i,9)  +
     >                biomass(i,10) +
     >                biomass(i,11) +
     >                biomass(i,12) +
     >                biomass(i,13) +
     >                biomass(i,14)  
c
c initial single-sided sai for upper and lower canopies
c
        sai(i,1)  =  0.050 * totlail(i)
        sai(i,2)  =  0.250 * totlaiu(i)
c
c fractional cover
c
c       fu(i) = wood / woodnorm
c
        fu(i) = (1.0 - exp(-wood)) / (1.0 - exp(-woodnorm))
c
        fl(i) = totlail(i) / 1.0
c
        fu(i) = max (0.25, min (0.975, fu(i)))
        fl(i) = max (0.25, min (0.975, fl(i)))
c
c initial lai for canopy physics
c
        lai(i,1) = totlail(i) / fl(i)
        lai(i,2) = totlaiu(i) / fu(i)
c
c specify canopy height parameters
c calculated as a function of only the vegetative fraction
c of each grid cell
c
c TET changed 6/26/02
c        zbot(i,1) =  0.05
c        ztop(i,1) =  max (0.25, lai(i,1) * 0.25)
        zbot(i,1) =  0.01
        ztop(i,1) =  0.50
c
c        zbot(i,2) =  ztop(i,1) + 1.0 
c        ztop(i,2) =  max (zbot(i,2) + 1.00, 
c     >                    2.50 * totbiou(i) / fu(i) * 0.75)
        zbot(i,2) =  2.0 
c CJK 1-16-03 as per TET        ztop(i,2) =  5.0 
        ztop(i,2) =  20.0 
c
  100 continue
c
c ************************************************************************
c assign some physical properties of vegetation
c ************************************************************************
c
c leaf optical properties were taken from Sellers et al., 1996
c and Bonan, 1995
c
*      rhoveg(1,1) = 0.10     ! vis leaf reflectance, lower story
*      rhoveg(1,2) = 0.10     ! vis leaf reflectance, upper story 
c
*      rhoveg(2,1) = 0.60     ! nir leaf reflectance, lower story
*      rhoveg(2,2) = 0.40     ! nir leaf reflectance, upper story
c
*      tauveg(1,1) = 0.07     ! vis leaf transmittance, lower story
*      tauveg(1,2) = 0.05     ! vis leaf transmittance, upper story
c
*      tauveg(2,1) = 0.25     ! nir leaf transmittance, lower story
*      tauveg(2,2) = 0.20     ! nir leaf transmittance, upper story
c
      if (iwheat .gt. 0) then
        chiflz = 0.65 
      else
        chiflz = -0.5        ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
      endif
c CJK chiflz =  0.0          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
      chifuz =  0.0          ! leaf orientation factors (-1 vertical, 0 random, 1 horizontal)
c
      oriev(1) = max (-chiflz, 0.)
      oriev(2) = max (-chifuz, 0.)
c
      orieh(1) = max ( chiflz, 0.)
      orieh(2) = max ( chifuz, 0.)
c
*      dleaf(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
*      dstem(1) = 0.10        ! linear dimensions for aerodynamic flux parameterization
c
*      dleaf(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
*      dstem(2) = 0.10        ! linear dimensions for aerodynamic flux parameterization
c
*      chu = ch2o *  2.0      ! heat capacity of upper leaves
*      chl = ch2o *  2.0      ! heat capacity of lower leaves
*      chs = ch2o * 50.0      ! heat capacity of stems
c
*      alaimu = 8.0           ! normalization constant for upper canopy aerodynamics
*      alaiml = 8.0           ! normalization constant for lower canopy aerodynamics
c
*      cleaf  = 0.01          ! constant in leaf-air aero transfer parameterization
*      cgrass = 0.01          ! constant in leaf-air aero transfer parameterization
*      cstem  = 0.01          ! constant in leaf-air aero transfer parameterization
c
*      wliqumax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
*      wliqsmax = 0.40        ! intercepted water capacity (mm h2o per unit leaf area)
*      wliqlmax = 0.20        ! intercepted water capacity (mm h2o per unit leaf area)
c
*      wsnoumax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
*      wsnosmax = 4.00        ! intercepted snow capacity (mm h2o per unit leaf area)
*      wsnolmax = 2.00        ! intercepted snow capacity (mm h2o per unit leaf area)
c
*      tdripu =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
*      tdrips =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
*      tdripl =  2.0 * 3600.0 ! decay time for intercepted liquid dripoff (sec)
c
*      tblowu = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
*      tblows = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
*      tblowl = 12.0 * 3600.0 ! decay time for snow blowoff (sec)
c
c ************************************************************************
c define rooting profiles
c ************************************************************************
c
c define rooting profiles based upon data published in:
c
c Jackson et al., 1996:  A global analysis of root distributions
c for terrestrial biomes, Oecologia, 108, 389-411.
c
c and
c
c Jackson et al., 1997:  A global budget for fine root biomass, 
c surface area, and nutrient contents, Proceedings of the National
c Academy of Sciences, 94, 7362-7366.
c
c rooting profiles are defined by the "beta" parameter
c
c beta1 is assigned to the lower vegetation layer (grasses and shrubs)
c beta2 is assigned to the upper vegetation layer (trees)
c
c according to Jackson et al. (1996, 1997), the values of beta
c typically fall in the following range
c
c note that the 1997 paper specifically discusses the distribution
c of *fine roots* (instead of total root biomass), which may be more
c important for water and nutrient uptake
c
c --------------                 ------------   ------------
c forest systems                 beta2 (1996)   beta2 (1997)
c --------------                 ------------   ------------
c tropical evergreen forest:        0.962          0.972
c tropical deciduous forest:        0.961          0.982
c temperate conifer forest:         0.976          0.980
c temperate broadleaf forest:       0.966          0.967
c all tropical/temperate forest:    0.970  
c boreal forest:                    0.943          0.943
c all trees:                                       0.976
c
c -------------------------      ------------   ------------
c grassland / shrub systems      beta1 (1996)   beta1 (1997)
c -------------------------      ------------   ------------
c tropical grassland / savanna:     0.972          0.972
c temperate grassland:              0.943          0.943
c all grasses:                      0.952          0.952
c schlerophyllous shrubs:           0.964          0.950
c all shrubs:                       0.978          0.975
c crops:                            0.961
c desert:                           0.975          0.970
c tundra:                           0.914
c
c --------------                 ------------
c all ecosystems                 beta  (1996)
c --------------                 ------------
c all ecosystems:                   0.966
c
c for global simulations, we typically assign the following
c values to the beta parameters
c
c beta1 = 0.950, which is typical for tropical/temperate grasslands
c beta2 = 0.970, which is typical for tropical/temperate forests
c
c however, these values could be (and should be) further refined
c when using the model for specific regions
c 
*       beta1 = 0.950  ! for lower layer herbaceous plants
*       beta2 = 0.975  ! for upper layer trees
c
c calculate total depth in centimeters
c
      totdepth = 0.0
c
      do 300 k = 1, nsoilay
        totdepth = totdepth + hsoi(k) * 100.0
 300  continue
c
c normalization factors
c
      frootnorm1 = 1. - beta1 ** totdepth
      frootnorm2 = 1. - beta2 ** totdepth
c
c calculate rooting profiles
c
      do 400 k = 1, nsoilay
c
        if (k.eq.1) then
c
          depth(k) = hsoi(k) * 100.0
c
          froot(k,1) = 1. - beta1 ** depth(k)
          froot(k,2) = 1. - beta2 ** depth(k)
c
        else
c
          depth(k) = depth(k-1) + hsoi(k) * 100.0
c
          froot(k,1) = (1. - beta1 ** depth(k)) - 
     >                 (1. - beta1 ** depth(k-1)) 
c
          froot(k,2) = (1. - beta2 ** depth(k)) - 
     >                 (1. - beta2 ** depth(k-1)) 
c
        endif
c
        froot(k,1) = froot(k,1) / frootnorm1
        froot(k,2) = froot(k,2) / frootnorm2
c
 400  continue
c
c return to main program
c
      return
      end
c
