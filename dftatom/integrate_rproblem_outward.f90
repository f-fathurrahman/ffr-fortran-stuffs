! Integrates the radial Schroedinger/Dirac equations outward for the given E
!
! Input parameters:
! l .... the quantum number l
! E .... the energy at which to integration the equation
! R(:) .... radial grid
! V(:) .... potential on the radial grid
! Z .... the nucleus charge
! c .... speed of light
! relat ... the relativistic mode (1, 2 or 3), see solve_radial_eigenproblem()
!           for more information
!
! Output parameters:
! P(:), Q(:) ..... The P and Q components.
! For Schroedinger equation P(r) = r*R(r), Q(r) = P'(r), for Dirac equation,
! P(r) = r*g(r), Q(r) = r*f(r), where "g" and "f" are the large and small
! components of the Dirac equation
!
! The components are not normalized.

!------------------------------------------------
SUBROUTINE integrate_rproblem_outward( &
   Nr, l, E, R, Rp, V, Z, c, relat, P, Q, imax )
!-----------------------------------------------
  IMPLICIT NONE 
  INTEGER :: Nr
  INTEGER, INTENT(in) :: l, relat, Z
  REAL(8), INTENT(in) :: E, c, R(Nr), Rp(Nr), V(Nr)
  REAL(8), INTENT(out) :: P(Nr), Q(Nr)
  INTEGER, INTENT(out) :: imax

  ! Suppress warning
  !INTEGER :: kappa
  REAL(8) :: ctmp
  ctmp = 1.d0*c

  IF( relat == 0 ) THEN 
    CALL schroed_outward_adams( Nr, l, Z, E, R, Rp, V, P, Q, imax )
  !
!  ELSEIF( relat == 1 ) THEN 
!    CALL stop_error("Scalar relativistic case not implemented")
!  !
!  ELSEIF( relat == 2 .or. relat == 3 ) THEN 
!    !
!    IF( relat == 3 ) THEN 
!      IF( l == 0 ) THEN 
!        CALL stop_error("for l=0 only spin up (relat=2) is allowed")
!      ENDIF 
!      kappa = l
!    !
!    ELSE 
!        kappa = -l-1
!    ENDIF
!    CALL dirac_outward_adams( c, kappa, Z, E, R, Rp, V, P, Q, imax )
  ELSE 
    !
    CALL stop_error("wrong value of relat or not yet supported")
  ENDIF 

END SUBROUTINE 

