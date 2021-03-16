! Integrates the following set of equations outwards:
! dy1/dx =           y2
! dy2/dx = C1 + C2 * y2
! The above system can be equivalently written as a second order ODE:
! y1'' - C2 * y1' = C1
! The coefficients C1 and C2 are passed in as arrays, so they can have any
! dependence on R. For example for a Poisson equation of the form:
! V'' + 2 * V' / R = -4*pi*n
! we would get C1 = -4*pi*n, C2 = -2/R
SUBROUTINE rk4_integrate3( Nr, R, y0, C1, C2, C1mid, C2mid, max_val, y1, y2, imax)
  !
  IMPLICIT NONE 
  INTEGER :: Nr
  REAL(8), INTENT(in) :: R(Nr) ! Grid
  REAL(8), INTENT(in) :: y0(2) ! Initial condition
  ! Coefficients C1 and C2 at grid points and midpoints:
  REAL(8), INTENT(in) :: C1(Nr), C2(Nr), C1mid(Nr-1), C2mid(Nr-1)
  ! Maximum value (if y1 > max_val, the integration stops)
  REAL(8), INTENT(in) :: max_val
  ! Solution y1 and y2
  REAL(8), INTENT(out) :: y1(Nr), y2(Nr)
  ! The integration stops at R(imax)
  INTEGER, INTENT(out) :: imax
  
  INTEGER :: i
  REAL(8) :: dym(2), dyt(2), yt(2), dydx(2), y(2)
  REAL(8) :: h

  y = y0
  y1(1) = y(1)
  y2(1) = y(2)
  DO i = 2,Nr
    !
    h = R(i) - R(i-1)
    !
    dydx(1) = y(2)
    dydx(2) = C1(i-1) + C2(i-1) * y(2)
    !
    yt = y + h/2 * dydx
    !
    dyt(1) = yt(2)
    dyt(2) = C1mid(i-1) + C2mid(i-1) * yt(2)
    !
    yt = y + h/2 * dyt
    !
    dym(1) = yt(2)
    dym(2) = C1mid(i-1) + C2mid(i-1) * yt(2)
    !
    yt = y + h * dym
    !
    dym = dyt + dym
    dyt(1) = yt(2)
    dyt(2) = C1(i) + C2(i) * yt(2)
    !
    y = y + h/6 * (dydx + dyt + 2*dym)
    !
    y1(i) = y(1)
    y2(i) = y(2)
    IF( abs(y(1)) >= max_val ) THEN 
      imax = i
      RETURN 
    ENDIF 
  ENDDO 
  imax = Nr
  RETURN 
END SUBROUTINE

