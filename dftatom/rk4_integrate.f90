!
! Integrates the following set of equations outwards:
! dy1/dx =                y2
! dy2/dx = C1 * y1 + C2 * y2
!
!--------------------------------------------------------------------------------
SUBROUTINE rk4_integrate( N, R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax )
!--------------------------------------------------------------------------------
  IMPLICIT NONE 
  INTEGER, INTENT(in) :: N
  REAL(8), INTENT(in) :: R(N) ! Grid
  REAL(8), INTENT(in) :: y0(2) ! Initial condition
  ! Coefficients C1 and C2 at grid points and midpoints:
  REAL(8), INTENT(in) :: C1(N), C2(N), C1mid(N-1), C2mid(N-1)
  ! Maximum value (if y1 > max_val, the integration stops)
  REAL(8), INTENT(in) :: max_val
  ! Solution y1 and y2
  REAL(8), INTENT(out) :: y1(N), y2(N)
  ! The integration stops at R(imax)
  INTEGER, INTENT(out) :: imax

  INTEGER :: i
  REAL(8) :: dym(2), dyt(2), yt(2), dydx(2), y(2)
  REAL(8) :: h

  y = y0
  y1(1) = y(1)
  y2(1) = y(2)
  DO i = 2, n
    ! step size
    h = R(i) - R(i-1)
    !
    dydx(1) = y(2)
    dydx(2) = C1(i-1) * y(1) + C2(i-1) * y(2)
    !
    yt = y + h/2 * dydx
    !
    dyt(1) = yt(2)
    dyt(2) = C1mid(i-1) * yt(1) + C2mid(i-1) * yt(2)
    !
    yt = y + h/2 * dyt
    !
    dym(1) = yt(2)
    dym(2) = C1mid(i-1) * yt(1) + C2mid(i-1) * yt(2)
    !
    yt = y + h * dym
    dym = dyt + dym
    !
    dyt(1) = yt(2)
    dyt(2) = C1(i) * yt(1) + C2(i) * yt(2)
    !
    y = y + h/6 * (dydx + dyt + 2*dym)

    y1(i) = y(1)
    y2(i) = y(2)
    IF( abs(y(1)) >= max_val ) THEN 
      imax = i
      RETURN 
    ENDIF 
  ENDDO 
  imax = n
  RETURN 
END SUBROUTINE 

