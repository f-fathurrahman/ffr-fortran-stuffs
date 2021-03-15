PROGRAM test_get_n_orb
  IMPLICIT NONE 
  INTEGER :: n_orb

  CALL get_n_orb( 11, n_orb )
  WRITE(*,*) 'n_orb = ', n_orb
END PROGRAM 
