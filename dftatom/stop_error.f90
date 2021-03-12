SUBROUTINE stop_error(msg)
  IMPLICIT NONE 
  CHARACTER(len=*) :: msg
  WRITE(*,*) msg
  STOP 1
END SUBROUTINE 
