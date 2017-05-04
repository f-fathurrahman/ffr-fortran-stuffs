MODULE m_simpleOps

CONTAINS 

SUBROUTINE increment(a, b)
  IMPLICIT NONE 
  INTEGER, INTENT(inout) :: a(:)
  INTEGER, INTENT(in) :: b
  INTEGER :: i, N

  N = size(a)
  DO i = 1,N
    a(i) = a(i) + b
  ENDDO 
END SUBROUTINE 

END MODULE 


!------------------------
PROGRAM increment_testCPU
!------------------------
  
  USE m_simpleOps
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N = 256
  INTEGER :: a(N), b

  a(:) = 1
  b = 3

  CALL increment(a,b)
  IF( any(a /= 4) ) THEN 
    WRITE(*,*) '!!!! Program Failed !!!!!'
  ELSE 
    WRITE(*,*) 'Program success!'
  ENDIF 

END PROGRAM 

