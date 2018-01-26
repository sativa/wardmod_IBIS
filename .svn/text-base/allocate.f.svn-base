c Bill Sacks
c 01.29.10

c ---------------------------------------------------------------------
      subroutine allocate_vars(niter)
c ---------------------------------------------------------------------
c
c Allocate variables
c
c Note: Dynamically allocated variables is a Fortran-90 feature. As
c such, this could break a Fortran-77-only compiler.
c
      use comatm, only : comatm_init
      implicit none
c
c Arguments
c
      integer niter  ! number of time iterations per day
      
      call comatm_init(niter)

      return
      end
