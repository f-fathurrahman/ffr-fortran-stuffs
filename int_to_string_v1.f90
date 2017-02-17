PROGRAM int_to_string
  IMPLICIT NONE
  INTEGER :: a, b
  CHARACTER(32) :: str

  a = 4
  b = 10

  WRITE(str,'(A,I1)') 'hello_', a
  WRITE(*,*) 'str = ', trim(str)

  WRITE(str,'(A,I2)') 'hello_', b
  WRITE(*,*) 'str = ', trim(str)

END PROGRAM
