!--------------------------------------------------
SUBROUTINE adams_interp_inward_implicit(N, y, i, r)
!--------------------------------------------------
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8), INTENT(in) :: y(N)
  INTEGER, INTENT(in) :: i
  REAL(8) :: r
  r = -(19*y(i) - 5*y(i+1) + y(i+2)) / 24d0
  RETURN 
END SUBROUTINE 

