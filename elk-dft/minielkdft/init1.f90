! Generates the $k$-point set and then allocates and initialises global
! variables which depend on the $k$-point set.
SUBROUTINE init1()
  USE modmain, ONLY: timeinit
  IMPLICIT NONE 
  ! local variables
  REAL(8) :: ts0, ts1

  CALL timesec(ts0)

  CALL init_kpoints()

  CALL init_Gk_vectors()

  CALL init_APW_LO()
 
  CALL init_dftu_k()

  CALL init_eigensystems()

  CALL timesec(ts1)
  timeinit = timeinit+ts1-ts0

  WRITE(*,*) '*** Leaving init1'

  RETURN 
END SUBROUTINE 

