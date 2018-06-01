PROGRAM test
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N = 5
  COMPLEX(8) :: z1(N), z2(N)
  COMPLEX(8) :: res, zdotc

  z1(:) = cmplx(1.d0,2.d0)
  z2(:) = cmplx(2.d0,3.d0)
  
  res = zdotc( N, z1, 1, z2, 1 )

  WRITE(*,*) 'res = ', res
END PROGRAM 

