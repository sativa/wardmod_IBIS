c
c ---------------------------------------------------------------------
      module comsno
c ---------------------------------------------------------------------
c
      use comgrid
      use compar, only : nsnolay

      implicit none
      save

      real 
     >  z0sno,               ! roughness length of snow surface (m)
     >  rhos,                ! density of snow (kg m-3)
     >  consno,              ! thermal conductivity of snow (W m-1 K-1)
     >  hsnotop,             ! thickness of top snow layer (m)
     >  hsnomin,             ! minimum total thickness of snow (m)
     >  fimin,               ! minimum fractional snow cover
     >  fimax                ! maximum fractional snow cover
c
      real 
     >  fi(npoi)             ! fractional snow cover
c
      real 
     >  tsno(npoi,nsnolay),  ! temperature of snow layers (K)
     >  hsno(npoi,nsnolay)   ! thickness of snow layers (m)
c
      end module comsno
