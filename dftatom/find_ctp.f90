! Finds the classical turning point for the potential 'V' and energy 'E'
! Classical turning point 'ctp' is defined as E = V(ctp)
! The function returns the integer index into the array V.
!--------------------------------
SUBROUTINE find_ctp(N, V, E, ctp)
!--------------------------------
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8), INTENT(in) :: V(N), E
  INTEGER :: i
  INTEGER :: ctp
  DO i = N, 1, -1
    IF( V(i) - E <= 0.d0 ) THEN 
      ctp = i
      RETURN 
    ENDIF 
  ENDDO 
  ctp = 0
  RETURN 
END SUBROUTINE 

