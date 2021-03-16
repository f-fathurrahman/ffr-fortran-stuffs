SUBROUTINE atom_lda( &
    Z, r_min, r_max, a, Nr, n_orb, no, lo, fo, ks_energies, E_tot, &
    R, Rp, V_tot, density, orbitals, reigen_eps, reigen_max_iter, &
    mixing_eps, mixing_alpha, &
    mixing_max_iter, perturb)

  USE m_dft_data, ONLY: dft_data_t
  USE m_states, ONLY: get_atomic_states_nonrel
  IMPLICIT NONE 
  !
  INTEGER, INTENT(in) :: Z
  REAL(8), INTENT(in) :: r_min, r_max, a
  INTEGER, INTENT(in) :: Nr
  INTEGER :: n_orb
  INTEGER, INTENT(out), TARGET :: no(n_orb), lo(n_orb)
  REAL(8), INTENT(out), TARGET :: fo(n_orb)
  REAL(8), INTENT(out), TARGET :: ks_energies(n_orb)
  REAL(8), INTENT(out) :: E_tot
  REAL(8), INTENT(out), TARGET :: R(Nr), Rp(Nr)
  REAL(8), INTENT(out), TARGET :: V_tot(Nr)
  REAL(8), INTENT(out), TARGET :: density(Nr)
  REAL(8), INTENT(out), TARGET :: orbitals(Nr, n_orb)
  REAL(8), INTENT(in) :: reigen_eps, mixing_eps, mixing_alpha
  INTEGER, INTENT(in) :: mixing_max_iter, reigen_max_iter
  LOGICAL, INTENT(in) :: perturb

  INTEGER, POINTER :: no_a(:), lo_a(:)
  REAL(8), POINTER :: fo_a(:)
  INTEGER :: i
  REAL(8), TARGET :: Emin_init(n_orb), Emax_init(n_orb) 
  REAL(8), TARGET :: V_h(Nr), V_xc(Nr), e_xc(Nr), V_coulomb(Nr), tmp(Nr)
  !
  TYPE(dft_data_t) :: d

  CALL get_atomic_states_nonrel( Z, no_a, lo_a, fo_a )
  no = no_a
  lo = lo_a
  fo = fo_a
  DEALLOCATE( no_a, lo_a, fo_a )

  CALL gen_rmesh_exp(r_min, r_max, a, Nr, R)
  CALL gen_drmesh_exp(r_min, r_max, a, Nr, Rp)

  !V_tot = -Z / R
  CALL get_hydrogen_energies( Z, n_orb, no, ks_energies )
  CALL calc_TF_potential( Nr, R, Z, .true., V_tot )
  !ks_energies = get_tf_energies(Z, no, fo)

  V_coulomb = -Z/R

  ! We allow a few unbounded states
  Emax_init = 10
  
  ! Use Hydrogen Schroedinger energies for the lower limit:
  DO i = 1,n_orb
    CALL calc_E_nl( 0.d0, no(i), lo(i), Z, 0, Emin_init(i) )
  ENDDO 
  ! For robustness, decrease Emin by 10%:
  Emin_init = 1.1d0 * Emin_init

  WRITE(*,*) 'Emin_init = ', Emin_init

  d%Z = Z
  d%R => R
  d%Rp => Rp
  d%rho => density
  d%V_h => V_h
  d%V_coulomb => V_coulomb
  d%V_xc => V_xc
  d%e_xc => e_xc
  d%V_tot => V_tot
  d%orbitals => orbitals
  d%alpha = mixing_alpha
  d%dirac = .false.
  d%perturb = perturb
  d%reigen_eps = reigen_eps
  d%reigen_max_iter = reigen_max_iter
  d%no => no
  d%lo => lo
  d%fo => fo
  d%ks_energies => ks_energies
  d%Emax_init => Emax_init
  d%Emin_init => Emin_init
  
  ! Start from an initial density instead:
  !d%rho = 1 / cosh(d%R)**2
  !d%rho = d%rho * d%Z / integrate(d%Rp, 4*pi*d%rho*d%R**2)
  !call rho2V(d)
  
  CALL scf_potmix_anderson( Nr, n_orb, d%V_tot - d%V_coulomb, mixing_max_iter, &
    .true., d, mixing_alpha, mixing_eps, tmp )
  
  E_tot = d%Etot
  ! Prints the energies:
  print *, "Ekin: ", d%Ekin
  print *, "Ecoul:", d%Ecoul
  print *, "Eenuc:", d%Eenuc
  print *, "Exc:  ", d%Exc
  print *, "E_tot:", d%Ekin + d%Ecoul + d%Eenuc + d%Exc
END SUBROUTINE 

