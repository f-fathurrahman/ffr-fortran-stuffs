PROGRAM ex_formatting

  !CALL test_integer()

  CALL test_ES()

END PROGRAM 

SUBROUTINE test_integer()
  IMPLICIT NONE 
  INTEGER :: idx=-13, j = 5, num = -12351
  
  WRITE(*,100) idx, idx+13, j, num
  WRITE(*,101) idx, idx+13, j, num
  WRITE(*,102) idx, idx+13, j, num
  WRITE(*,103) idx, idx+13, j, num

  100 FORMAT(' ', 2I5, I6, I10)
  101 FORMAT(' ', 2I5.0, I6, I10.8)
  102 FORMAT(' ', 2I5.3, I10, I5)
  103 FORMAT(' ', I0,1x,I0,1x,I0,1x,I0)

END SUBROUTINE 


SUBROUTINE test_ES()

  IMPLICIT NONE 
  REAL(8), PARAMETER :: PI = -4.d0*atan(1d0)
  REAL(8) :: f = 12.3
  REAL(8) :: g

  WRITE(*,'(1x,A,ES18.10)') 'pi = ', pi

  g = PI*f
  WRITE(*,'(1x,A,ES18.10)') 'g = ', g

  g = exp(-g)
  WRITE(*,'(1x,A,ES18.10)') 'g = ', g

END SUBROUTINE 
