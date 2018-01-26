c ---------------------------------------------------------------------
      module combgc
c ---------------------------------------------------------------------
c
c ----------------------------------------
c Soil biogeochemistry parameters for ibis
c ----------------------------------------
c
      use comgrid

      implicit none
      save

c constants for calculating fraction of litterall in structural
c metabolic and lignified (resistant) fractions
c
      real 
     >     lig_frac,   ! split of lignified litter material between protected/non-protected slow OM pools
     >     fbsom,      ! protected biomass as a fraction of total soil organic C from Verberne et al., 1990
     >     effac,      ! efficiency of microbial biomass reincorporated into biomass pool. 
*                      ! (From NCSOIL parameterizations; Molina et al., 1983)
     >     cnr(10),    ! C:N ratios of substrate pools and biomass for leaves and roots.
*                      ! Values from Parton et al., 1987 and Whitmore and Parry, 1988
     >     fmxcpool,       ! maximum fraction allowed in resistant fraction (Verbene 1997)
     >     rconst,     ! constant defined as 1200 (from Verbene 1997 equations)
     >     h20         ! C leaching fudge factor

c
c decay constants for c pools
c
      real
     >          klm,          ! leaf metabolic litter 
     >          kls,          ! leaf structural litter
     >          kll,          ! leaf lignin
     >          krm,          ! root metabolic litter
     >          krs,          ! root structural litter
     >          krl,          ! root lignin
     >          kwm,          ! woody metabolic litter
     >          kws,          ! woody structural litter
     >          kwl,          ! wood  lignin
     >          kbn,          ! microbial biomass --> nonprotected om 
     >          kbp,          ! microbial biomass --> protected om
     >          knb,          ! nonprotected om   --> biomass
     >          kns,          ! nonprotected om   --> passive c 
     >          kpb,          ! protected om      --> biomass
     >          kps,          ! protected om      --> passive c
     >          ksb           ! passive c         --> biomass

c
c efficiencies for microbial decomposition
c
      real
     >          ylm,          ! leaf metabolic litter decomposition 
     >          yls,          ! leaf structural litter decomposition
     >          yll,          ! leaf lignin
     >          yrm,          ! root metabolic litter decomposition
     >          yrs,          ! root structural litter decomposition
     >          yrl,          ! root lignin
     >          ywm,          ! woody metabolic litter decomposition
     >          yws,          ! woody structural litter decomposition
     >          ywl,          ! wood lignin
     >          ybn,          ! microbial biomass to nonprotected om
     >          ybp,          ! microbial biomass to protected om
     >          ynb,          ! nonprotected om to biomass
     >          yns,          ! nonprotected om to passive  c
     >          ypb,          ! protected om to biomass
     >          yps,          ! protected om to passive c
     >          ysb           ! passive c to biomass
c
      end module combgc
