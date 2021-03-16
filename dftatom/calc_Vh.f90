SUBROUTINE calc_Vh( Nr, R, Rp, rho, V )
  IMPLICIT NONE 
  INTEGER :: Nr
  REAL(8), intent(in) :: R(Nr), Rp(Nr), rho(Nr)
  REAL(8) :: V(Nr)
  CALL rpoisson_outward_pc( Nr, R, Rp, rho, V )
END SUBROUTINE 

