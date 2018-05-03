! File myCG.f90
MODULE myCG
CONTAINS
SUBROUTINE myCGf(v, f, nz, ny, dx, dy, eps)
implicit none
integer, intent(in) :: nz, ny
real (kind=8), intent(in) :: dx,dy
real (kind=8), dimension(0:(nz+1),0:(ny+1)), intent(inout):: v
real (kind=8), dimension(0:(nz+1),0:(ny+1)), intent(in) :: f
real (kind=8), intent(in) :: eps
real (kind=8), dimension(0:(nz+1),0:(ny+1)) :: q, r, rho
real (kind=8), external :: dnrm2, ddot
real (kind=8) :: r2, beta, alpha, r2New
integer :: i,k
integer :: ntot, nout
! Tell f2py what is input and output
!f2py intent(in) :: f, nz, ny, dx, dy, eps
!f2py intent(out) :: V
nout=0
ntot = (nz+2)*(ny+2)

! CG stuff happens here
! v is updated and returned in the python environment
END SUBROUTINE

END MODULE
