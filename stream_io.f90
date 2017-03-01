PROGRAM writeUstream
  IMPLICIT NONE
  INTEGER :: myvalue = 12345, mypos
  OPEN(UNIT=11, FILE="ustream.demo", STATUS="replace", ACCESS="STREAM")
  WRITE(11) "first"
  WRITE(11) "second"
  INQUIRE(UNIT=11, POS=mypos)
  PRINT *, "Myvalue will be written at position ", mypos
  WRITE(11) myvalue
  CLOSE(UNIT=11)
END PROGRAM
