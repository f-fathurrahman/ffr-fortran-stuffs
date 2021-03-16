SUBROUTINE rho_to_V( d )
  USE m_dft_data, ONLY: dft_data_t
  ! Calculates V_xc, V_h and V_tot from rho.
  ! Assumes that d%rho is normalized.
  TYPE(dft_data_t), INTENT(inout) :: d
  INTEGER :: Nr, n_orb

  Nr = size(d%R)
  n_orb = size(d%fo)

  CALL calc_Vxc( Nr, d%R, d%rho, d%dirac, d%c, d%e_xc, d%V_xc )
  CALL calc_Vh( Nr, d%R, d%Rp, d%rho, d%V_h )
  !
  CALL calc_Etot( n_orb, d%fo, d%ks_energies, d%V_tot, d%V_h, d%V_coulomb, d%e_xc, &
    Nr, d%R, d%Rp, d%rho, d%Ekin, d%Ecoul, d%Eenuc, d%Exc, d%Etot )
  !
  d%V_tot = d%V_coulomb + d%V_xc + d%V_h

  RETURN 
END SUBROUTINE 

