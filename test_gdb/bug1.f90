PROGRAM bug

  IMPLICIT  NONE 
  REAL :: last
  REAL :: c(10)
  INTEGER :: p

  ! Initialise c with successive integer values.
  DO p = 1,10
    c(p) = p
  ENDDO 

  ! Calculate and print ratios of successive integers.
  last = 0.0
  DO p=1,10
    CALL divide(last, c(p))
    last = c(p)
  ENDDO

END PROGRAM 

SUBROUTINE divide(d,e)
  IMPLICIT NONE 
  REAL :: d,e
  WRITE(*,*) e/d
END SUBROUTINE 

