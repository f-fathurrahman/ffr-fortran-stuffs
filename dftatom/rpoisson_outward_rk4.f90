! Solves V''(r) + 2/r*V'(r) = -4*pi*density (in this form) with initial
! conditions V(R(1)) = V0 and V'(R(1)) = V0d using 4th order Runge-Kutta method
! by rewriting it to the equivalent first order ODEs (see the documentation
! of rk4_integrate3() for more details). Returns V (value) and Vd (derivative).
SUBROUTINE rpoisson_outward_rk4( Nr, density, density_mid, R, V0, V0d, V, Vd )
  IMPLICIT NONE 
  INTEGER :: Nr
  REAL(8), INTENT(in) :: density(Nr), density_mid(Nr-1) ! density at the
  ! mesh points R and at midpoints
  REAL(8), INTENT(in) :: R(Nr) ! radial grid
  REAL(8), INTENT(in) :: V0  ! value of V(r) at r=R(1)
  REAL(8), INTENT(in) :: V0d ! value of V'(r) at r=R(1)
  REAL(8), INTENT(out) :: V(Nr), Vd(Nr)  ! V(r) and V'(r) on the grid R
  !
  REAL(8), PARAMETER :: pi=4.d0*atan(1.d0)

  REAL(8) :: y(2)
  REAL(8) :: C1(Nr), C2(Nr)
  REAL(8) :: C1mid(Nr-1), C2mid(Nr-1)
  REAL(8) :: Rmid(Nr-1)
  INTEGER :: imax

  y(1) = V0
  y(2) = V0d

  C1 = -4.d0*pi*density
  C2 = -2.d0 / R
  Rmid = ( R(:size(R)-1) + R(2:) ) / 2
  C1mid = -4.d0*pi*density_mid
  C2mid = -2.d0 / Rmid
  CALL rk4_integrate3( Nr, R, y, C1, C2, C1mid, C2mid, 1.0d10, V, Vd, imax )

  IF( imax /= size(R) ) CALL stop_error("Poisson solver diverged.")
END SUBROUTINE 

