PROGRAM simple_loop
  
  IMPLICIT NONE 

  INTEGER :: i
  INTEGER, PARAMETER :: N=20
  REAL(8) :: A(N), B(N)

  A(:) = 1.d0
  B(:) = 2.1d0

!$OMP PARALLEL DO
  DO i=2,N
    WRITE(*,*) 'i = ', i
    B(i) = 0.5d0 * ( A(i) + A(i-1) )
  ENDDO 
!$OMP END PARALLEL DO

END PROGRAM 
