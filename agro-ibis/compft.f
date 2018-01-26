c ---------------------------------------------------------------------
      module compft
c     last modified by C. Kucharik 10.04.02
c ---------------------------------------------------------------------
c
c -------------------------
c plant parameters for ibis
c -------------------------

      use compar, only : npft

      implicit none
      save

      real*4 
     >       tau15,          ! co2/o2 specificity ratio at 15 degrees C (dimensionless)
     >       kc15,           ! co2 kinetic parameter (mol/mol)
     >       ko15,           ! o2 kinetic parameter (mol/mol) 
     >       cimax,          ! maximum value for ci (needed for model stability)
     >       alpha3,         ! intrinsic quantum efficiency for C3 plants (dimensionless)
     >       theta3,         ! photosynthesis coupling coefficient for C3 plants (dimensionless)
     >       thetac3,        ! photosynthesis coupling coefficient for C3 crops
     >       beta3,          ! photosynthesis coupling coefficient for C3 plants (dimensionless)
     >       alpha4,         ! intrinsic quantum efficiency for C4 plants (dimensionless)
     >       beta4,          ! photosynthesis coupling coefficient for C4 plants (dimensionless)
     >       theta4,         ! photosynthesis coupling coefficient for C4 plants (dimensionless)
     >       thetac4,        ! photosynthesis coupling coefficient for C4 crops
     >       betac3,         ! photosynthesis coupling coefficient for C3 crops 
     >       betac4,         ! photosynthesis coupling coefficient for C4 crops
     >       vmax_pft(npft)  ! nominal vmax of top leaf at 15 C (mol-co2/m**2/s) [not used]

c broadleaf
      real*4 
     >       gammaub,        ! leaf respiration coefficient
     >       coefmub,        ! 'm' coefficient for stomatal conductance relationship
     >       coefbub,        ! 'b' coefficient for stomatal conductance relationship
     >       gsubmin         ! absolute minimum stomatal conductance

c conifer
      real*4 
     >       gammauc,        ! leaf respiration coefficient
     >       coefmuc,        ! 'm' coefficient for stomatal conductance relationship  
     >       coefbuc,        ! 'b' coefficient for stomatal conductance relationship  
     >       gsucmin         ! absolute minimum stomatal conductance

c shrub
      real*4 
     >       gammals,        ! leaf respiration coefficient
     >       coefbls,        ! 'b' coefficient for stomatal conductance relationship 
     >       coefmls,        ! 'm' coefficient for stomatal conductance relationship 
     >       gslsmin         ! absolute minimum stomatal conductance

c c3 grass
      real*4 
     >       gammal3,        ! leaf respiration coefficient
     >       coefml3,        ! 'm' coefficient for stomatal conductance relationship 
     >       coefbl3,        ! 'b' coefficient for stomatal conductance relationship 
     >       gsl3min         ! absolute minimum stomatal conductance


c c4 grass
      real*4 
     >       gammal4,        ! leaf respiration coefficient
     >       coefml4,        ! 'm' coefficient for stomatal conductance relationship
     >       coefbl4,        ! 'b' coefficient for stomatal conductance relationship
     >       gsl4min         ! absolute minimum stomatal conductance

c c3 crop
      real*4 
     >       gammac3,        ! leaf respiration coefficient
     >       coefmc3,        ! 'm' coefficient for stomatal conductance relationship 
     >       coefbc3,        ! 'b' coefficient for stomatal conductance relationship 
     >       gsc3min         ! absolute minimum stomatal conductance

c c4 crop
      real*4 
     >       gammac4,        ! leaf respiration coefficient
     >       coefmc4,        ! 'm' coefficient for stomatal conductance relationship
     >       coefbc4,        ! 'b' coefficient for stomatal conductance relationship
     >       gsc4min         ! absolute minimum stomatal conductance

c crop
      real*4 
     >     lotemp(npft),      ! low temperature threshold in tempvm equation
     >     hitemp(npft),      ! high temperature threshold in tempvm equation 
     >     drought(npft),     ! crop sensitivity to drought parameter 
     >     f1(npft),          ! constant used in tempvm equations 
     >     f2(npft)           ! constant used in tempvm equations
     
c general pft     
      real*4 
     >       tauleaf(npft),   ! foliar biomass turnover time constant (years)
     >       tauroot(npft),   ! fine root biomass turnover time constant (years)
     >       tauwood(npft),   ! wood biomass turnover time constant (years)
     >       tauwood0(npft),  ! wood biomass turnover time constant (years) 
     >       plai_init(4,15), ! initial total LAI for each vegtype (used in iniveg)
     >       plaiupper,       ! Potental LAI of upper canopy (uniform initial vegetation)
     >       plailower,       ! Potental LAI of lower canopy (uniform initial vegetation)
     >       xminlai,         ! Minimum LAI for each existing PFT
     >       sapfrac_init,    ! Initial value of sapwood fraction used for all woody PFTs
     >       beta1,          ! parameter for Jackson rooting profile (lower canopy)
     >       beta2           ! parameter for Jackson rooting profile (upper canopy)

c more general pft
      real*4 
     >       TminL(npft),    ! Absolute minimum temperature -- lower limit (upper canopy PFTs)
     >       TminU(npft),    ! Absolute minimum temperature -- upper limit (upper canopy PFTs)
     >       Twarm(npft),    ! Temperature of warmest month (lower canopy PFTs)
     >       GDD(npft)       ! minimum GDD needed (base 5 C for upper canopy PFTs, 
                             ! base 0 C for lower canopy PFTs)
      
      end module compft


 

