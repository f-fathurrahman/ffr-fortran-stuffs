! from Brainerd
PROGRAM trig

  IMPLICIT NONE 

  REAL(8) :: x = 0.4d0

  x = 0.1d0 + this_image()/10.d0

  SELECT CASE( this_image() )

    CASE(1)
      WRITE(*,*) 'sin(',x,') = ', sin(x)

    CASE(2)
      WRITE(*,*) 'cos(',x,') = ', cos(x)

    CASE(3)
      WRITE(*,*) 'tan(',x,') = ', tan(x)

    CASE DEFAULT
      WRITE(*,*) 'Default case: not calculating anything'

  END SELECT 

END PROGRAM 

