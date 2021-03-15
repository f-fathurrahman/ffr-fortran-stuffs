MODULE m_dft_data

! Contains the 'dft_data_t' type used in the DFT routines.
! This data type stores mesh, potential, atomic configuration, orbitals
! and other parameters of the DFT problem.

IMPLICIT NONE 

TYPE dft_data_t
  REAL(8), DIMENSION(:), POINTER :: R, Rp, V_coulomb, V_h, V_xc, e_xc, V_tot, rho
  REAL(8), DIMENSION(:, :), POINTER :: orbitals
  REAL(8) :: reigen_eps, alpha, c
  INTEGER :: Z, scf_iter, reigen_max_iter
  ! If .true., we are solving the Dirac equation
  LOGICAL :: dirac
  ! Triples of (n, l, f), where "n", "l" are quantum numbers and "f" is the
  ! occupation number:
  INTEGER, DIMENSION(:), POINTER :: no, lo, so
  REAL(8), DIMENSION(:), POINTER :: fo
  REAL(8), POINTER :: ks_energies(:), Emax_init(:), Emin_init(:)
  ! Total energy (and its parts):
  REAL(8) :: Ekin, Ecoul, Exc, Eenuc, Etot
  LOGICAL :: perturb ! use perturbation acceleration in the eigenproblem
END TYPE 

END MODULE 

