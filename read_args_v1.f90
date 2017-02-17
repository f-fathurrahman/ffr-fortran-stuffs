PROGRAM read_args
  IMPLICIT NONE
  REAL(8) :: arg1
  INTEGER :: arg2
  CHARACTER(10) :: arg3
  CHARACTER(20) :: buffer

  CALL getarg(1,buffer)
  READ(buffer,*) arg1

  CALL getarg(2,buffer)
  READ(buffer,*) arg2

  CALL getarg(3,buffer)
  READ(buffer,*) arg3

  WRITE(*,*) 'First argument : ', arg1
  WRITE(*,*) 'Second argument: ', arg2
  WRITE(*,*) 'Third argument : ', arg3
END PROGRAM

