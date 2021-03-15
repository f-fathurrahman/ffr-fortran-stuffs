SUBROUTINE get_n_orb( Z, n_orb )
  USE m_states, ONLY: get_atomic_states_nonrel
  IMPLICIT NONE 
  INTEGER :: Z
  INTEGER :: n_orb
  
  INTEGER, POINTER, DIMENSION(:) :: no, lo
  REAL(8), POINTER, DIMENSION(:) :: fo

  CALL get_atomic_states_nonrel(Z, no, lo, fo)
  n_orb = size(no)
  DEALLOCATE(no, lo, fo)
END SUBROUTINE 

