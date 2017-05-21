! from Brainerd
PROGRAM trig

  IMPLICIT NONE 
  CHARACTER(*), PARAMETER :: STRFMT = '(1x,A,2F18.10)'


  REAL(8) :: x = 0.4d0

  x = 0.1d0 + this_image()/10.d0

  SELECT CASE( this_image() )

    CASE(1)
      WRITE(*,STRFMT) 'x, sin(x) = ', x, sin(x)

    CASE(2)
      WRITE(*,STRFMT) 'x, cos(x) = ', x, cos(x)

    CASE(3)
      WRITE(*,STRFMT) 'x, tan(x) = ', x, tan(x)

    CASE DEFAULT
      WRITE(*,*) 'Default case: not calculating anything'

  END SELECT 

END PROGRAM 

