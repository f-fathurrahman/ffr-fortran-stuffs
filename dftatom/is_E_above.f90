!---------------------------------------------
SUBROUTINE is_E_above(n, l, nods_actual, res)
!---------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n, l, nods_actual
  LOGICAL :: res 
  res = nods_actual > n-l-1
  RETURN 
END SUBROUTINE 
