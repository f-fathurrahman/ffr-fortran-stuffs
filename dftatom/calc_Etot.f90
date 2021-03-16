SUBROUTINE calc_Etot( n_orb, fo, ks_energies, V_in, V_h, V_coulomb, e_xc, &
        Nr, R, Rp, n, &
        T_s, E_ee, E_en, EE_xc, Etot)
  IMPLICIT NONE 
  INTEGER :: Nr, n_orb
  ! This is a variational, quadratically convergent form of total energy
  REAL(8), INTENT(in) :: fo(n_orb), ks_energies(n_orb) ! occupations, energies
  REAL(8), INTENT(in) :: V_in(Nr) ! Total input effective potential
  REAL(8), INTENT(in) :: V_h(Nr) ! Hartree energy, solution of Poiss. eq.
  REAL(8), INTENT(in) :: V_coulomb(Nr) ! Coulomb inter. -Z/r  (negative)
  REAL(8), INTENT(in) :: e_xc(Nr) ! XC density
  REAL(8), INTENT(in) :: R(Nr), Rp(Nr), n(Nr) ! Radial grid, number density (positive)
  REAL(8), INTENT(out) :: Etot ! Total energy
  REAL(8), INTENT(out) :: T_s, E_ee, E_en, EE_xc ! Parts of the total energy

  REAL(8), PARAMETER :: pi = 4.d0*atan(1.d0)
  REAL(8) :: rho(Nr)
  REAL(8) :: E_c, E_band !, Exc2
  REAL(8) :: ss

  rho = -n  ! why???

  E_band = sum(fo * ks_energies)
  CALL integrate_trapz_7( Nr, Rp, V_in * rho * R**2, ss )
  T_s = E_band + 4*pi * ss

  CALL integrate_trapz_7( Nr, Rp, V_h * rho * R**2, ss )
  E_ee = -2*pi * ss
  !
  CALL integrate_trapz_7( Nr, Rp, (-V_coulomb)*rho*R**2, ss )
  E_en =  4*pi * ss
  !
  E_c = E_ee + E_en

  CALL integrate_trapz_7( Nr, Rp, e_xc * rho * R**2, ss )
  EE_xc = -4*pi * ss

  Etot = T_s + E_c + EE_xc

  RETURN 
END SUBROUTINE 

