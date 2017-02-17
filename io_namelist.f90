PROGRAM io_namelist

  IMPLICIT NONE
  
  REAL :: i_real4
  REAL(8) :: i_real8
  INTEGER :: i_integer
  CHARACTER(10) :: i_char
  LOGICAL :: i_bool

  INTEGER, PARAMETER :: IUNIT=11

  NAMELIST /indata/ i_real4, i_real8, i_integer, i_char, i_bool

  OPEN(IUNIT, file="IN_io_namelist", STATUS="old")
  READ(IUNIT, nml=indata)
  CLOSE(IUNIT)

  WRITE(*, nml=indata)

  WRITE(*,*) 'Program ended'

END PROGRAM 
