PROGRAM test_atom_lda
  IMPLICIT NONE 
  
  ! Atomic number:
  INTEGER :: Z
  ! Mesh parameters:
  REAL(8), PARAMETER :: r_min = 1d-7, r_max = 50.d0, a = 2.7d6
  INTEGER :: Nr

  integer :: i, n_orb
  CHARACTER, PARAMETER :: l_names(0:3) = (/ "s", "p", "d", "f" /)
  REAL(8) :: E_tot
  REAL(8), PARAMETER :: reigen_eps = 1d-10
  REAL(8), PARAMETER :: mixing_eps = 5d-9
  INTEGER, ALLOCATABLE, dimension(:) :: no, lo
  REAL(8), ALLOCATABLE, dimension(:) :: fo, ks_energies
  REAL(8), ALLOCATABLE :: orbitals(:, :)
  REAL(8), ALLOCATABLE :: R(:), Rp(:), V_tot(:), density(:)

  Z = 5
  Nr = 5500 + 1
  CALL get_n_orb(Z, n_orb)

  ALLOCATE( ks_energies(n_orb) )
  ALLOCATE( no(n_orb) )
  ALLOCATE( lo(n_orb) )
  ALLOCATE( fo(n_orb) )

  ALLOCATE( orbitals(Nr, n_orb) )
  ALLOCATE( R(Nr) )
  ALLOCATE( V_tot(Nr) )
  ALLOCATE( density(Nr) )
  ALLOCATE( Rp(Nr) )
  
  CALL atom_lda( Z, r_min, r_max, a, Nr, n_orb, no, lo, fo, ks_energies, E_tot, &
    R, Rp, V_tot, density, orbitals, reigen_eps, 100, mixing_eps, &
    0.5d0, 200, .true. )

  WRITE(*,*) "Z=", Z, "Nr = ", Nr
  WRITE(*,'(1x,A,F18.10)') "E_tot=", E_tot
  WRITE(*,*) "state    E            occupancy"
  DO i = 1, size(ks_energies)
    WRITE(*,"(1x,I1, A, ' ', F18.6, '   ', F6.3)") no(i), l_names(lo(i)), &
          ks_energies(i), fo(i)
  ENDDO

  DEALLOCATE( ks_energies )
  DEALLOCATE( no )
  DEALLOCATE( lo )
  DEALLOCATE( fo )

  DEALLOCATE( orbitals )
  DEALLOCATE( R )
  DEALLOCATE( V_tot )
  DEALLOCATE( density )
  DEALLOCATE( Rp )

END PROGRAM 

