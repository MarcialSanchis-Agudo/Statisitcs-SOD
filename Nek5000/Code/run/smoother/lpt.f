#include "experimental/meshsmoother.f"

c-----------------------------------------------------------------------
C
C  USER SPECIFIED ROUTINES:
C
C     - boundary conditions
C     - initial conditions
C     - variable properties
C     - local acceleration for fluid (a)
C     - forcing function for passive scalar (q)
C     - general purpose routine for checking errors etc.
C
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)

      parameter (init_trans=0)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'


      ffx=0.0
      ffy=0.0
      ffz=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'


      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'


      ux=1.0
      uy=0.0
      uz=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux=1.0
      uy=0.0
      uz=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'
c
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

      call meshsmoother 

ccc   The default parameters for the smoother are:
c      parameter(nbc=1)         !number of boundary conditions
c      character*3 dcbc(nbc)
c      save        dcbc
c      data        dcbc /'W  '/  !BCs listed here

c      idftyp = 0      !distance function - 0 -> exponential, 1-> tanh
c      alpha = 15.     !Input for wall distance function 
c      beta  = 0.1     !

c      nouter = 50      !total loops around laplacian and optimizer smoothing
c      nlap = 20        !number of laplacian iterations in each loop
c      nopt = 20        !number of optimization iterations in each loop

c      mtyp = 1         !metric type
ccc
ccc    NOTE:The user can modify these parameters and run the smoother 
ccc    with modified parameters as:
ccc    call smoothmesh(mtyp,nouter,nlap,nopt,nbc,dcbc,idftyp,alpha,beta)

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'
      return
      end
c-----------------------------------------------------------------------

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end