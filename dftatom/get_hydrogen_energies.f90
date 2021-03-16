SUBROUTINE get_hydrogen_energies(Z, n_orb, no, E)
  IMPLICIT NONE 
  INTEGER :: n_orb
  INTEGER :: Z, no(n_orb)
  REAL(8) :: E(n_orb)
  INTEGER :: i
  DO i = 1,n_orb
    E(i) = -1.d0 * Z**2 / ( 2*no(i)**2 )
  ENDDO 
  RETURN 
END SUBROUTINE 
