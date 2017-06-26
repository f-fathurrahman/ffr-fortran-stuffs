SUBROUTINE sub1()
  IMPLICIT NONE 
  INTEGER :: i
  DO i = 1,1000
    WRITE(*,*) 'I am called from test_flush: ', i
  ENDDO 
END SUBROUTINE 

SUBROUTINE sub2()
  IMPLICIT NONE 
  INTEGER :: i
  DO i = 1,1000
    WRITE(*,*) 'I am called again from test_flush: ', i
  ENDDO 
END SUBROUTINE 

PROGRAM test_flush
  IMPLICIT NONE 

  CALL sub1()
  CALL flush(6)
  CALL sub2()

END PROGRAM 
