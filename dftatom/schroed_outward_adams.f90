! Integrates the Schrodinger equation outward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form outwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = C * P
! where C = 2*(V-E) + l*(l+1)/R**2
!
! Returns P(r), Q(r).
!
! It tranforms the problem to a uniform mesh 1 <= t <= N + 1, defined by
! R(t) and Rp(t) = dR/dt and Rpp(t) = d^2R/dt^2:
!
! u(t)   = P(R(t))
! up(t)  = du/dt  = dP/dR * dR/dt             = Q * Rp
! upp(t) = dup/dt = dQ/dR * Rp^2 + Q * dRp/dt = C*P*Rp^2 + Q * Rpp
!
! So
!
! up  = u * Rp
! upp = u * C * Rp^2 + up * Rpp / Rp
!
! For example if Rp = al*R, Rpp = al^2 * R, we get:
!
! up  = u * al * R
! upp = u * C * al^2 * R^2 + up * al

!-------------------------------------------------------------------
SUBROUTINE schroed_outward_adams( N, l, Z, E, R, Rp, V, P, Q, imax )
!-------------------------------------------------------------------
  IMPLICIT NONE 
  INTEGER, INTENT(in) :: N
  INTEGER, INTENT(in) :: l
  INTEGER, INTENT(in) :: Z
  REAL(8), INTENT(in) :: E
  REAL(8), INTENT(in) :: R(N), Rp(N)
  REAL(8), INTENT(in) :: V(N)
  REAL(8), INTENT(out) :: P(N), Q(N)
  REAL(8), PARAMETER :: max_val = 1.0d6
  INTEGER, INTENT(out) :: imax ! The integration was carried to R(imax)
  
  REAL(8), DIMENSION(N) :: C, u1, u2, u1p, u2p, Vmid  ! XXX automatic arrays
  INTEGER :: i
  REAL(8) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp
  REAL(8) :: xx1, xx2

  IF( N < 5 ) THEN
    CALL stop_error("size(R) < 5")
  ENDIF 
  
  IF( .not. (size(R) == size(Rp) .and. size(R) == size(V) .and. &
      size(R) == size(P) .and. size(R) == size(P) .and. size(R) == size(Q))) then
    CALL stop_error("Array sizes mismatch")
  ENDIF

  WRITE(*,*) 'N = ', N
  WRITE(*,*) 'size(R) = ', size(R)
  
  C = ( l*(l+1)/R**2 + 2 * (V-E) )
  
  ! Get initial values from RK4
  !Vmid(:3) = get_midpoints(R(:4), V(:4))
  CALL get_midpoints( 4, R(:4), V(:4), Vmid(:3) )
  CALL integrate_rschroed_rk4( 4, l, Z, E, R(:4), V(:4), Vmid(:3), u1(:4), u2(:4), imax )

  !WRITE(*,*) 'imax after rk4 = ', imax
  !STOP 'after rk4'

  !u1(1:4) = R(1:4) ** (l+1)
  !u2(1:4) = (l+1) * R(1:4) ** l
  u1p(1:4) = Rp(1:4)          * u2(1:4)
  u2p(1:4) = Rp(1:4) * C(1:4) * u1(1:4)

  DO i = 4,N-1
    u1p(i) = Rp(i)        * u2(i)
    u2p(i) = Rp(i) * C(i) * u1(i)
    !
    CALL adams_interp_outward_implicit(N, u1p, i, xx1)
    CALL adams_interp_outward_implicit(N, u2p, i, xx2)
    u1_tmp  = u1(i) + xx1
    u2_tmp  = u2(i) + xx2

    lambda = 9.d0/24d0
    Delta = 1 - lambda**2 * C(i+1) * Rp(i+1)**2
    M(1, 1) = 1 / Delta
    M(2, 1) = lambda * C(i+1) * Rp(i+1) / Delta
    M(1, 2) = lambda * Rp(i+1) / Delta
    M(2, 2) = 1 / Delta

    u1(i+1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i+1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
    IF( abs(u1(i+1)) >= max_val .or. abs(u2(i+1)) >= max_val ) THEN 
      P = u1
      Q = u2
      imax = i
      RETURN 
    ENDIF 
  ENDDO 
  !
  P = u1
  Q = u2
  imax = size(R)
  !
  RETURN 
END SUBROUTINE 

