SUBROUTINE read_input( filename )

  IMPLICIT NONE 
  CHARACTER(*) :: filename
  INTEGER, PARAMETER :: IUNIT = 2909
  CHARACTER(64) :: line(3)
  INTEGER :: i
  INTEGER :: ios
  
  WRITE(*,*) 'filename = ', filename

  OPEN(unit=IUNIT, file=filename, status='old', action='read', iostat=ios)
  IF( ios /= 0 ) THEN 
    WRITE(*,*) 'ERROR reading file: ', trim(filename)
    STOP 
  ENDIF 

  DO WHILE( .TRUE. )

    READ(IUNIT,*,END=1988) line

    DO i = 1,3
      WRITE(*,*) i, trim(line(i))
    ENDDO 

  ENDDO 

  1988 CONTINUE 
  CLOSE(IUNIT)

END SUBROUTINE 

 ! split a string into 2 either side of a delimiter token
SUBROUTINE split_string(instring, string1, string2, delimiter)
  IMPLICIT NONE 
    CHARACTER(128) :: instring
    CHARACTER(1) :: delimiter
    CHARACTER(128), INTENT(OUT):: string1,string2
    INTEGER :: idx

    instring = TRIM(instring)

    idx = SCAN(instring,delimiter)
    WRITE(*,*) 'idx = ', idx
    string1 = instring(1:idx-1)
    WRITE(*,*) instring(1:idx-1)
    string2 = instring(idx+1:)
    WRITE(*,*) instring(idx+1:)

END SUBROUTINE


