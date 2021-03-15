!---------------------
PROGRAM test_gen_rmesh
!---------------------
  IMPLICIT NONE 
  REAL(8) :: r_min, r_max, a
  INTEGER :: Nr
  REAL(8), ALLOCATABLE :: rmesh(:)
  INTEGER :: i
  
  Nr = 12
  r_min = 0.d0
  r_max = 50.d0
  a = 1d9
  !a = 1.d0
  ALLOCATE( rmesh(Nr) )
  
  CALL gen_rmesh_exp(r_min, r_max, a, Nr, rmesh)

  DO i = 1,Nr
    WRITE(*,'(1x,I8,ES18.10)') i, rmesh(i)
  ENDDO

  DEALLOCATE(rmesh)
END PROGRAM 

