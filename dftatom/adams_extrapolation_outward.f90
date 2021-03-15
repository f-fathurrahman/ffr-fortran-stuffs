!-----------------------------------------------------
SUBROUTINE adams_extrapolation_outward( N, y, i, res )
!-----------------------------------------------------
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: y(N)
  INTEGER :: i
  REAL(8) :: res
  res = +( 55*y(i) - 59*y(i-1) + 37*y(i-2) - 9*y(i-3) ) / 24
END SUBROUTINE 

