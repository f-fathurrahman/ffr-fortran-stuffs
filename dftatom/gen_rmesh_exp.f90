!--------------------------------------------------
SUBROUTINE gen_rmesh_exp(r_min, r_max, a, Nr, rmesh)
!--------------------------------------------------
  IMPLICIT NONE 
  REAL(8), INTENT(in) :: r_min, r_max
  REAL(8), INTENT(in) :: a
  INTEGER, INTENT(in) :: Nr
  REAL(8) :: rmesh(Nr)

  INTEGER :: i, N
  REAL(8) :: alpha, beta

  N = Nr - 1

  IF ( a < 0.d0 ) THEN 
    CALL stop_error("mesh_exp: a > 0 required")
  !
  ! a == 1  
  !
  ELSEIF( abs(a - 1.d0) < tiny(1.d0) ) THEN 
    alpha = (r_max - r_min) / N
    DO i = 1,Nr
      rmesh(i) = alpha*(i - 1.0d0) + r_min
    ENDDO 
  !
  ! normal case
  ELSE 
    IF( N > 1 ) THEN 
        beta = log(a)/(N-1)
        alpha = (r_max - r_min) / (exp(beta*N) - 1)
        DO i = 1, N+1
          rmesh(i) = alpha * ( exp(beta*(i-1)) - 1 ) + r_min
        ENDDO 
    ELSEIF( N == 1 ) THEN 
        rmesh(1) = r_min
        rmesh(2) = r_max
    ELSE 
      CALL stop_error("gen_rmesh_exp: N >= 1 required")
    ENDIF 
  ENDIF 
END SUBROUTINE 
