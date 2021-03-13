! Integrates the Schrodinger equation outward using Runge-Kutta 4th order method
!
! It integrates the Schroedinger equation in the R(r) form outwards:
! R'' = -2/R * R' + (2*(V-E) + l*(l+1)/R**2) * R
!
! Returns P(r), Q(r), where these are defined as:
! P(r) = r*R(r)
! Q(r) = P'(r) = r * R'(r) + R(r)
! where R(r) is the radial solution.

!--------------------------------------------------------------------
SUBROUTINE integrate_rschroed_rk4(N, l, Z, E, R, V, Vmid, P, Q, imax)
!--------------------------------------------------------------------

  IMPLICIT NONE 
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: l
  INTEGER, INTENT(in) :: Z
  REAL(8), INTENT(in) :: E
  REAL(8), INTENT(in) :: R(N)
  REAL(8), INTENT(in) :: V(N), Vmid(N-1)
  REAL(8), INTENT(out) :: P(N), Q(N)
  INTEGER, INTENT(out) :: imax ! The integration was carried to R(imax)
  
  REAL(8), PARAMETER :: max_val = 1.d6
  REAL(8) :: y0(2)
  REAL(8) :: y1(N), y2(N)  ! XXX automatic arrays
  REAL(8), DIMENSION(N) :: C1, C2
  REAL(8), DIMENSION(N-1) :: C1mid, C2mid
  REAL(8), DIMENSION(N-1) :: Rmid
  
  ! y(:) are values of the components (2) of the equation, in our case:
  ! y(1) = R, y(2) = R'
  ! dydx(:) are the derivatives of the components (2) of the equation

  IF( l == 0 ) THEN 
    y0(1) = 1-Z*R(1)
    y0(2) = -Z
  ELSE 
    y0(1) = R(1)**l
    y0(2) = l*R(1)**(l-1)
  ENDIF 

  IF( size(V) /= size(Vmid) + 1 ) THEN  
    CALL stop_error("Vmid size is wrong")
  ENDIF 

  ! Evaluate array values needed for RK4
  C1 = 2*(V-E) + l*(l+1)/R**2
  C2 = -2/R
  Rmid = (R(:size(R)-1) + R(2:)) / 2
  C1mid = 2*(Vmid-E) + l*(l+1)/Rmid**2
  C2mid = -2/Rmid

  CALL rk4_integrate(N, R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)

  P(:imax) = y1(:imax)*R(:imax) ! P(r) = r * R(r)
  Q(:imax) = y2(:imax)*R(:imax) + y1(:imax) ! Q(r) = P'(r) = r * R'(r) + R(r)

END SUBROUTINE 
