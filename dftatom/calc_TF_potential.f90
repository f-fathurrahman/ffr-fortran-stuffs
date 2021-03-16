! Generalized Thomas-Fermi atomic potential
SUBROUTINE calc_TF_potential( Nr, R, Z, cut, V )
  IMPLICIT NONE 
  INTEGER :: Nr
  REAL(8), INTENT(in) :: R(Nr)
  INTEGER, INTENT(in) :: Z
  LOGICAL, INTENT(in) :: cut ! Cut the potential
  REAL(8) :: x(Nr), Z_eff(Nr), V(Nr)
  REAL(8) :: alpha, beta, gamma1
  REAL(8), PARAMETER :: pi=4.d0*atan(1.d0)

  x = R * (128*Z/(9*pi**2)) ** (1.d0/3)
  ! Z_eff(x) = Z * phi(x), where phi(x) satisfies the Thomas-Fermi equation:
  !   phi'' = phi**(3/2) / sqrt(x)
  ! with boundary conditions:
  !   phi(0)  = 1
  !   phi(oo) = 0
  ! There is no analytic solution, but one can solve this approximately. We use:
  ! http://arxiv.org/abs/physics/0511017
  alpha = 0.7280642371d0
  beta = -0.5430794693d0
  gamma1 = 0.3612163121d0
  !
  Z_eff = Z * ( 1 + alpha*sqrt(x) + &
    beta*x*exp(-gamma1*sqrt(x)))**2 * exp(-2*alpha*sqrt(x))

  ! This keeps all the eigenvalues of the radial problem negative:
  IF (.not. cut ) WHERE (Z_eff < 1.d0) Z_eff = 1.d0
  
  V = -Z_eff / r

  RETURN 
END SUBROUTINE 

