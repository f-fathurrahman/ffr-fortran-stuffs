! Solves the equation V''(r) + 2/r*V'(r) = -4*pi*rho
!
! Uses predictor corrector Adams method.
!
! It rewrites it to the equivalent system of first order ODEs on a uniform
! grid:
!   u1 = V
!   u2 = V'
!   u1p = u2 * Rp
!   u2p = -(4*pi*rho + 2*u2/r) * Rp
! and integrates outward using Adams method. The initial conditions are:
!   V (R(1)) = u1(1) = 4*pi * \int r * rho(r) dr
!   V'(R(1)) = u2(1) = 0
SUBROUTINE rpoisson_outward_pc(Nr, R, Rp, rho, V)
  IMPLICIT NONE 
  INTEGER :: Nr
  REAL(8), INTENT(in) :: R(Nr), Rp(Nr), rho(Nr)
  REAL(8) :: V(Nr)
  !
  REAL(8), PARAMETER :: pi=4.d0*atan(1.d0)
  !
  REAL(8), ALLOCATABLE :: u1(:), u2(:), u1p(:), u2p(:)
  INTEGER :: i, it
  INTEGER, PARAMETER :: max_it = 2
  REAL(8) :: rho_mid(3), xx1, xx2, xx

  ALLOCATE( u1(Nr) )
  ALLOCATE( u2(Nr) )
  ALLOCATE( u1p(Nr) )
  ALLOCATE( u2p(Nr) )

  CALL get_midpoints( 4, R(:4), rho(:4), rho_mid )  
  
  CALL integrate_trapz_7( Nr, Rp, rho*R, xx )
  CALL rpoisson_outward_rk4( 4, rho(:4), rho_mid, R(:4), &
    4*pi*xx, 0.d0, u1(:4), u2(:4) )

  u1p(:4) = u2(:4) * Rp(:4)
  u2p(:4) = -(4*pi*rho(:4) + 2*u2(:4)/R(:4)) * Rp(:4)

  DO i = 4, Nr-1
    CALL adams_extrapolation_outward( Nr, u1p, i, xx1 )
    u1(i+1) = u1(i) + xx1
    CALL adams_extrapolation_outward( Nr, u2p, i, xx2 )
    u2(i+1) = u2(i) + xx2
    !
    DO it = 1, max_it
      u1p(i+1) = +Rp(i+1) * u2(i+1)
      u2p(i+1) = -Rp(i+1) * ( 4*pi*rho(i+1) + 2*u2(i+1)/R(i+1) )
      CALL adams_interp_outward( Nr, u1p, i, xx1 )
      u1(i+1)  = u1(i) + xx1
      CALL adams_interp_outward( Nr, u2p, i, xx2 )
      u2(i+1)  = u2(i) + xx2
    ENDDO 
  ENDDO 
  V = u1

  DEALLOCATE( u1 )
  DEALLOCATE( u2 )
  DEALLOCATE( u1p )
  DEALLOCATE( u2p )

END SUBROUTINE 

