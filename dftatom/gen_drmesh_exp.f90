!----------------------------------------------------
SUBROUTINE gen_drmesh_exp(r_min, r_max, a, N, drmesh)
!----------------------------------------------------
  IMPLICIT NONE 
  REAL(8), INTENT(in) :: r_min
  REAL(8), INTENT(in) :: r_max
  REAL(8), INTENT(in) :: a
  INTEGER, INTENT(in) :: N
  REAL(8) :: drmesh(N+1)
  
  INTEGER :: i
  REAL(8) :: alpha, beta
  
  IF( a < 0.d0 ) THEN 
    CALL stop_error("gen_drmesh_exp: a > 0 required")
  !
  ELSEIF( abs(a - 1) < tiny(1.d0) ) THEN 
    CALL stop_error("gen_drmesh_exp: a == 1 not implemented")
    ! FIXME: use ones ??
  !
  ELSE 
    IF( N > 1 ) THEN 
      beta = log(a)/(N-1)
      alpha = (r_max - r_min) / (exp(beta*N) - 1)
      DO i = 1, N+1
        drmesh(i) = alpha * beta * exp(beta*(i-1))
      ENDDO 
    ELSE 
      CALL stop_error("mesh_exp_deriv: N > 1 required")
    ENDIF 
  !
  ENDIF 
END SUBROUTINE 

