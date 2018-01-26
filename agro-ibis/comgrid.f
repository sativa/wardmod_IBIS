c ---------------------------------------------------------------------
      module comgrid
c     parameters defining the ibis grid
c ---------------------------------------------------------------------
c WJS (02.11.09): I separated out these grid-related parameters from compar.h
c because they are likely to change from user to user and from run to run.
c
      integer nlon,    ! longitude dimension of domain
     >        nlat,    ! latitude dimension of domain
     >        npoi     ! total number of land points
c
      real xres,       ! longitude resolution (in degrees)
     >     yres        ! latitude resolution (in degrees)
c
c --------------------------
c typical ibis configuration
c --------------------------
c
c
c global 10.0 x 10.0 TEST grid: 36 by 18 array, with 147 land points:
c
c     parameter (nlon = 36,
c    >           nlat = 18, 
c    >           npoi = 147,
c    >           xres = 10.00,
c    >           yres = 10.00)
c
c global 0.5 x 0.5  grid:  720 by 360 array, with 58920 land points:
c
c      parameter (nlon = 720,
c     >           nlat = 360, 
c     >           npoi = 58920,
c     >           xres = 0.50,
c     >           yres = 0.50)
c
c
c global 1.0 x 1.0  grid:  360 by 180 array, with 14545 land points:
c
c      parameter (nlon = 360,
c     >           nlat = 180, 
c     >           npoi = 14545,
c     >           xres = 1.00, 
c     >           yres = 1.00)
c
c global 2.0 x 2.0  grid:  180 by  90 array, with  3728 land points:
c
c      parameter (nlon = 180,
c     >           nlat = 90, 
c     >           npoi = 3728,
c     >           xres = 2.00,
c     >           yres = 2.00)
c
c global 4.0 x 4.0  grid:   90 by  45 array, with   912 land points:
c
c     parameter (nlon = 90,
c    >           nlat = 45, 
c    >           npoi = 912,
c    >           xres = 4.00,
c    >           yres = 4.00)
c
c U.S. only 0.5 x 0.5 grid:   140 by 56 array, with   4838 land points:
c
c      parameter (nlon = 140,
c     >           nlat = 56,
c    >           npoi = 1707,
c     >           npoi = 3295, 
c    >           npoi = 4, 
c     >           xres = 0.50,
c     >           yres = 0.50)
c
c zedx-based 5min grid: 960 by 336 array, with 115562 (unextrap) or 159508 (extrap) land points:

c      parameter (nlon = 960,
c     >           nlat = 336,
c     >           npoi = 702,
cc     >           npoi = 115562,
c     >           xres = 0.08333333333,
c     >           yres = 0.08333333333)


c user-defined grid mode
c
c
      parameter (nlon = 960,
     >           nlat = 336,
     >           npoi = 1,  !1088,
     >           xres = 0.08333333333,
     >           yres = 0.08333333333)
c
      end module comgrid
