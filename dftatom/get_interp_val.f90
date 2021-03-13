! returns the interpolated value of f(x), where f is defined on the
! grid R
! input parameters
! R ... grid
! n ... length of the grid
! f ... function defined on the grid R
! x ... point at which to interpolate
! i ... the integer number of the radial point at which "x" lies
! output parameters
! V ... the interpolated function value
SUBROUTINE get_interp_val(f, x, R, n, V, i)
  IMPLICIT NONE 
  INTEGER, INTENT(in) :: n
  REAL(8), INTENT(in) :: f(n), x, R(n)
  REAL(8), INTENT(out) :: V
  INTEGER, INTENT(in) :: i
  
  INTEGER :: j1, j2, n1, n2
  
  IF( n < 4 ) CALL stop_error("get_interp_val: n >= 4 required")
  
  j1 = i - 1
  j2 = j1 + 1
  
  n1 = j1 - 1
  n2 = j2 + 1

  IF( n1 < 1 ) THEN 
      n2 = n2 - n1 + 1
      n1 = 1
  ENDIF 

  IF( n2 > n ) THEN 
      n1 = n1 - ( n2 - n)
      n2 = n
  ENDIF 
  
  CALL calc_interp( x, r(n1:n2), n2-n1+1, f(n1:n2), V )
END SUBROUTINE 

