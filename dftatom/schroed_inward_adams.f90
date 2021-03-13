! Integrates the Schrodinger equation inward using Adams 4th order method
!
! It integrates the Schroedinger equation in the P(r), Q(r) form inwards using
! predictor-corrector:
! P' = Q
! Q' = P'' = (2*(V-E) + l*(l+1)/R**2) * P
!
! Returns P(r), Q(r).
!
! Important note: This only works for exponential meshes, where the mesh
! marameter (as defined by mesh_exp()) is equal to "a = (rmax/rmin)**((N-1)/N)".
! Otherwise it will produce wrong answer.
SUBROUTINE schroed_inward_adams(Nr, l, E, R, Rp, V, P, Q, imin)
  IMPLICIT NONE 
  INTEGER :: Nr
  INTEGER, INTENT(in) :: l
  REAL(8), INTENT(in) :: E
  REAL(8), INTENT(in) :: R(Nr), Rp(Nr)
  REAL(8), INTENT(in) :: V(Nr)
  REAL(8), INTENT(out) :: P(Nr), Q(Nr)
  REAL(8), PARAMETER :: max_val = 1.0d300
  INTEGER, INTENT(out) :: imin

  REAL(8), ALLOCATABLE :: C(:), u1(:), u2(:), u1p(:), u2p(:) ! automatic arrays
  INTEGER :: i, i_max
  REAL(8) :: lambda, Delta, M(2, 2), u1_tmp, u2_tmp
  REAL(8) :: R_max
  REAL(8) :: xx1, xx2

  WRITE(*,*) 'in schroed_inward_adams: Nr = ', Nr
  WRITE(*,*) 'size(R) = ', size(R)

  allocate( C(Nr) )
  allocate( u1(Nr) )
  allocate( u2(Nr) )
  allocate( u1p(Nr) )
  allocate( u2p(Nr) )

  C = (l*(l+1)/R**2 + 2 * (V-E))

  i_max = Nr - 4
  IF( i_max < 4 ) CALL stop_error("imax < 4")
  
  lambda = sqrt(-2*E)
  ! We require that exp(-lambda*(R-R(1)) ~ epsilon(1.0_dp),
  ! if we start further from the origin, it might sometimes blow up.
  ! if we start closer, we might not get as precise asymptotic.
  ! It follows that R ~ R(1)-log(epsilon(1.0_dp)) / lambda
  R_max = R(1) - log( epsilon(1.d0) )/lambda
  DO WHILE( R(i_max) > R_max )
    IF( i_max == 2 ) THEN 
      CALL stop_error("Can't start the inward integration")
    ENDIF 
    i_max = i_max - 1
  ENDDO 

  write(*,*) 'Pass here 57 in schroed_inward_adams'

  u1(i_max:i_max+4) = exp(-lambda * (R(i_max:i_max+4)-R(1)))
  u2(i_max:i_max+4) = - lambda * u1(i_max:i_max+4)
  u1p(i_max:i_max+4) = Rp(i_max:i_max+4)                * u2(i_max:i_max+4)
  u2p(i_max:i_max+4) = Rp(i_max:i_max+4) * C(i_max:i_max+4) * u1(i_max:i_max+4)

  write(*,*) 'Pass here 64 in schroed_inward_adams'

  DO i = i_max,2,-1
    u1p(i) = Rp(i)        * u2(i)
    u2p(i) = Rp(i) * C(i) * u1(i)
    CALL adams_interp_inward_implicit( Nr, u1p, i, xx1 )
    u1_tmp  = u1(i) + xx1
    CALL adams_interp_inward_implicit( Nr, u2p, i, xx2 )
    u2_tmp  = u2(i) + xx2

    lambda = -9.d0 / 24
    Delta = 1 - lambda**2 * C(i-1) * Rp(i-1)**2
    M(1, 1) = 1 / Delta
    M(2, 1) = lambda * C(i-1) * Rp(i-1) / Delta
    M(1, 2) = lambda * Rp(i-1) / Delta
    M(2, 2) = 1 / Delta

    u1(i-1) = M(1, 1) * u1_tmp + M(1, 2) * u2_tmp
    u2(i-1) = M(2, 1) * u1_tmp + M(2, 2) * u2_tmp
    IF( abs(u1(i-1)) >= max_val .or. abs(u2(i-1)) >= max_val ) THEN 
      P = u1
      Q = u2
      P(i_max+1:) = 0
      Q(i_max+1:) = 0
      imin = i
      RETURN 
    ENDIF 
  ENDDO 
  P = u1
  Q = u2
  P(i_max+4:) = 0
  Q(i_max+4:) = 0
  imin = 1

  deallocate( C )
  deallocate( u1 )
  deallocate( u2 )
  deallocate( u1p )
  deallocate( u2p )

  RETURN 
END SUBROUTINE 
