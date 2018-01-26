c ---------------------------------------------------------------------
      module comdiag
c     last update 11/18/01 CJK 
c ---------------------------------------------------------------------

c
c This holds the parameters and global variables for diagnostics
c
      implicit none
      save
c
      integer nvars, nfiles
      parameter (nvars=128, nfiles=10)   ! nvars = # diagnostic variables
c
      integer ldiag(nvars,nfiles)      ! chosen diagnostic variables
c
c
      integer diagstart(nfiles),      ! year diagnostic output begins
     >        diagend(nfiles),        ! year diagnostic output stops
     >        ndiagpt(nfiles),        ! point in an npoi array
     >        nfreq(nfiles)           ! frequency of diagnostic output
c
      end module comdiag
