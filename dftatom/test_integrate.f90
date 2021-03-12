SUBROUTINE init_func(N, r, f)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: r(N)
  REAL(8) :: f(N)
  INTEGER :: i

  DO i = 1,N
    f(i) = exp(-0.2d0*r(i))
  ENDDO 
  RETURN 
END SUBROUTINE 


PROGRAM test_integrate
  IMPLICIT NONE 
  REAL(8) :: r_min, r_max, a
  INTEGER :: N
  REAL(8), ALLOCATABLE :: rmesh(:), drmesh(:), f(:)
  REAL(8) :: s
  
  N = 500
  r_min = 0.d0
  r_max = 200.d0
  a = 1d9
  !a = 1.d0
  ALLOCATE(rmesh(N+1))
  ALLOCATE(drmesh(N+1))
  ALLOCATE(f(N+1))

  CALL gen_rmesh_exp(r_min, r_max, a, N, rmesh)
  CALL gen_drmesh_exp(r_min, r_max, a, N, drmesh)

  CALL init_func(N+1, rmesh, f)

  WRITE(*,*) 'calculating integral ...'
  CALL integrate_trapz_7(N+1, drmesh, f, s)
  WRITE(*,*) 's = ', s

  DEALLOCATE(f)
  DEALLOCATE(drmesh)
  DEALLOCATE(rmesh)

  WRITE(*,*) 'Pass here'

END PROGRAM 
