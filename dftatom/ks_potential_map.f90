SUBROUTINE ks_potential_map( Nr, V, i, d, F )
  USE m_dft_data, ONLY: dft_data_t
  IMPLICIT NONE 
  ! Calculates the residual R(V) = KS(V) - V
  INTEGER :: Nr
  REAL(8), INTENT(in) :: V(Nr)
  INTEGER, INTENT(in) :: i
  TYPE(dft_data_t), INTENT(inout) :: d
  REAL(8) :: F(Nr)

  d%scf_iter = i
  ! We converge upon the non-Coulombic part of the potential:
  d%V_tot = d%V_coulomb + V
  CALL V_to_rho( d )
  CALL rho_to_V( d )
  F = d%V_xc + d%V_h - V
  RETURN 
END SUBROUTINE 

