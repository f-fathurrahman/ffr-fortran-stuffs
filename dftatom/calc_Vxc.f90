!----------------------------------------------------
SUBROUTINE calc_Vxc( Nr, R, rho, relat, c, exc, Vxc )
!----------------------------------------------------
  IMPLICIT NONE 
  INTEGER :: Nr
  REAL(8), INTENT(in) :: R(Nr)
  REAL(8), INTENT(in) :: rho(Nr) ! charge density
  LOGICAL, INTENT(in) :: relat ! .true. return RLDA, otherwise LDA
  REAL(8), INTENT(in) :: c 
  REAL(8), INTENT(out) :: Vxc(Nr), exc(Nr)

  INTEGER :: i
  DO i = 1,Nr
    CALL calc_Vxc_scalar( rho(i), relat, c, exc(i), Vxc(i) )
  ENDDO 
  RETURN 
END SUBROUTINE 

!--------------------------------------------------------
SUBROUTINE calc_Vxc_scalar( n, relat, c_light, exc, Vxc )
!--------------------------------------------------------
  IMPLICIT NONE 
  ! Calculates XC LDA density and potential from the charge density "n".
  REAL(8), INTENT(in) :: n ! charge density (scalar)
  REAL(8), INTENT(in) :: c_light ! speed of light
  LOGICAL, INTENT(in) :: relat ! if .true. returns RLDA, otherwise LDA
  REAL(8), INTENT(out) :: exc ! XC density
  REAL(8), INTENT(out) :: Vxc ! XC potential

  REAL(8), PARAMETER :: y0 = -0.10498d0
  REAL(8), PARAMETER :: b = 3.72744d0
  REAL(8), PARAMETER :: c = 12.9352d0
  REAL(8), PARAMETER :: A = 0.0621814d0
  REAL(8), PARAMETER :: pi=4.d0*atan(1.d0)

  REAL(8) :: Q, rs, y, ec, ex, Vc, Vx, beta, mu, R, S

  IF( n < epsilon(1.d0) ) THEN 
    exc = 0.d0
    Vxc = 0.d0
    RETURN 
  ENDIF 

  Q = sqrt(4*c - b**2)
  rs = (3/(4*pi*n))**(1.d0/3)
  y = sqrt(rs)
  ec = A/2 * (log(y**2/get_Y(y, b, c)) + 2*b/Q * atan(Q/(2*y+b))  &
   - b*y0/get_Y(y0, b, c) * ( &
            log((y-y0)**2 / get_Y(y, b, c)) &
            + 2*(b+2*y0) / Q * atan(Q/(2*y+b)) &
          ) )
  Vc = ec - A/6 * (c*(y-y0)-b*y0*y)/((y-y0)*get_Y(y, b, c))
  ex = -3/(4*pi) * (3*pi**2*n)**(1.d0/3)
  Vx = 4*ex/3

  IF( relat ) THEN 
    beta = -4*pi*ex / (3*c_light)
    mu = sqrt(1 + beta**2)
    R = 1 - 3*((beta * mu - log(beta + mu)) / (beta**2))**2 / 2
    S = 3 * log(beta + mu) / (2*beta * mu) - 1.d0/2
    ex = ex * R
    Vx = Vx * S
  ENDIF 
  exc = ex + ec
  Vxc = Vx + Vc

CONTAINS 

REAL(8) FUNCTION get_Y( y, b, c )
  REAL(8), INTENT(in) :: y, b, c
  get_Y = y**2 + b*y + c
END FUNCTION 

END SUBROUTINE 

