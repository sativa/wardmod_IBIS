c
c ---------------------------------------------------------------------
      module comwork
c     last update 2/02/10 Bill Sacks
c ---------------------------------------------------------------------
c
c this module holds work space arrays and the longitude and latitude
c vectors.  NOTE: The value of ndim3 should be nlon*nlat times the
c largest 3rd dimension in compar (max of nband,nsoilay,nsnolay,npft).
c The value ndim4 should be the largest of
c (nband,nsoilay,nsnolay,npft,nlon,nlat).
c
      use comgrid
      use compar

      implicit none
      save

      integer ndim2, ndim3, ndim4, max3d
      real OCEAN
c
      parameter ( OCEAN=9.e+20 ) ! the value to use for non-land points
      parameter ( ndim2 = nlon*nlat, 
     >            ndim4 = max(nlon, nlat, nband, nsoilay, nsnolay, npft), 
     >            max3d = max(nband, nsoilay, nsnolay, npft, 12),  ! WJS 02.02.10: largest 3rd dimension in compar; 12 included in the max expression for the sake of monthly arrays
     >  ndim3 = nlon*nlat*max3d ) 
c
      character aname*21       ! use to store short names
c
      integer lonindex (npoi), ! i index of nth point in land+sea array (i,j)
     >        latindex (npoi)  ! j index of nth point in land+sea array (i,j)
c
      real lonscale (nlon),    ! longitude of nth point in degrees east
     >     latscale (nlat)     ! latitude of nth point in degrees morth
c
      real work (ndim2)        ! work space big enough for one full grid
      real work2d (nlon, nlat) ! 2-d work space big enough for one full grid
c
c 

      real cdummy (ndim3)      ! work space big enough for npft grids
c
      end module comwork
