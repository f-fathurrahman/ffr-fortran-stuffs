!----------------------------------------------
SUBROUTINE adams_interp_outward( N, y, i, res )
!----------------------------------------------
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: y(N)
  INTEGER :: i
  REAL(8) :: res
  res = +( 9*y(i+1) + 19*y(i) - 5*y(i-1) + y(i-2) ) / 24
  RETURN 
END SUBROUTINE 

