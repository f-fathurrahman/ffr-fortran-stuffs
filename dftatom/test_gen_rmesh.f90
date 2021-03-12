!---------------------
PROGRAM test_gen_rmesh
!---------------------
  IMPLICIT NONE 
  REAL(8) :: r_min, r_max, a
  INTEGER :: N
  REAL(8), ALLOCATABLE :: rmesh(:)
  INTEGER :: i
  
  N = 11
  r_min = 0.d0
  r_max = 50.d0
  a = 1d9
  !a = 1.d0
  ALLOCATE( rmesh(N+1) )
  
  CALL gen_rmesh_exp(r_min, r_max, a, N, rmesh)

  DO i = 1,N+1
    WRITE(*,'(1x,I8,ES18.10)') i, rmesh(i)
  ENDDO

  DEALLOCATE(rmesh)
END PROGRAM 

