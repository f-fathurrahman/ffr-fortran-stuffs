SUBROUTINE calc_interp(t, X, n, Y, val)
  IMPLICIT NONE 
  ! Interpolates the x, y points to estimate the value at the point "t".
  ! val (out) ...... the estimated function value at the point "t".
  ! t (in).......... the "x" value at which to estimate the value
  ! X (array, in) .. the x-values of the points
  ! Y (array, in) .. the y-values of the points
  ! n (in) ......... the length of the arrays X and Y
  REAL(8), INTENT(in) :: t
  INTEGER, INTENT(in) :: n
  REAL(8), INTENT(in) :: X(n), Y(n)
  REAL(8), INTENT(out) :: val
  REAL(8) :: f, denum
  INTEGER :: i, j
  val = 0
  DO j = 1, n
    f = 1
    denum = 1
    DO i=1, n
      IF (i == j) CYCLE 
      f = f * ( t - X(i) )
      denum = denum * ( X(j) - X(i) )
    ENDDO 
    val = val + Y(j) * f / denum
  ENDDO
END SUBROUTINE 

