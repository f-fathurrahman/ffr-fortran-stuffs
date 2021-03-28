!---------------------------------------------------------------
SUBROUTINE integrate_rproblem_inward(Nr, l, E, R, Rp, V, c, &
    relat, P, Q, imin)
!---------------------------------------------------------------
  IMPLICIT NONE 
  INTEGER :: Nr
  INTEGER :: l, relat
  REAL(8) :: E, c, R(Nr), Rp(Nr), V(Nr)
  ! Outputs
  REAL(8) :: P(Nr), Q(Nr)
  INTEGER :: imin
  
  !INTEGER :: kappa
  REAL(8) :: ctmp
  ctmp = c

  WRITE(*,*) 'DEBUG: integrate_rproblem_inward, Nr = ', Nr

  IF( relat == 0 ) THEN 
    CALL schroed_inward_adams(Nr, l, E, R, Rp, V, P, Q, imin)
!  else if (relat == 1) then
!    call stop_error("Scalar relativistic case not implemented")
!  else if (relat == 2 .or. relat == 3) then
!    if (relat == 3) then
!      if (l == 0) then
!        call stop_error("for l=0 only spin up (relat=2) is allowed")
!      end if
!      kappa = l
!    else
!      kappa = -l-1
!    end if
!    call dirac_inward_adams(c, kappa, E, R, Rp, V, P, Q, imin)
  ELSE
    CALL stop_error("Wrong value of relat.")
  ENDIF 
END SUBROUTINE 
