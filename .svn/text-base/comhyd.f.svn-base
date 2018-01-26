c
c ---------------------------------------------------------------------
      module comhyd
c ---------------------------------------------------------------------
c
      use comgrid
      use compar,only: nsoilay

      implicit none
      save

      real 
     >  ginvap(npoi),                ! total evaporation rate from all intercepted h2o (kg_h2o m-2 s-1)
     >  gsuvap(npoi),                ! total evaporation rate from surface (snow/soil) (kg_h2o m-2 s-1)
     >  gtrans(npoi),                ! total transpiration rate from all vegetation canopies (kg_h2o m-2 s-1)
     >  gtransu(npoi),               ! transpiration from upper canopy (kg_h2o m-2 s-1)
     >  gtransl(npoi),               ! transpiration from lower canopy (kg_h2o m-2 s-1)
     >  grunof(npoi),                ! surface runoff rate (kg_h2o m-2 s-1)
     >  gdrain(npoi),                ! drainage rate out of bottom of <isoilay> soil layer (kg_h2o m-2 s-1)
     >  gdrain_layer(npoi,nsoilay),  ! drainage rate out of bottom of each soil layer   (kg_h2o m-2 s-1)
     >  gadjust(npoi),               ! h2o flux due to adjustments in subroutine wadjust (kg_h2o m-2 s-1)
     >  wtot(npoi)                   ! total amount of water stored in snow, soil, puddles, and on vegetation (kg_h2o)
c
      end module comhyd

