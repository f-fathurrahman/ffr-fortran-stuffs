PROGRAM add

  IMPLICIT NONE 
  INTEGER, PARAMETER :: DP=8
  REAL(8) :: x(1000)[*], y(1000)[*], z(1000)[*]
  INTEGER :: i

  IF( this_image() == 1 ) THEN 
    WRITE(*,*) 'Program running with ', num_images(), ' images'
  ENDIF 

  x(:) = 1.0_DP
  y(:) = 1.342_DP

  z(:) = x(:) + y(:)

  IF( this_image() == 1 ) THEN 
    DO i = 1,5
      WRITE(*,'(1x,I4,F20.15)') i, z(i)
    ENDDO 
    WRITE(*,*) 'Pass here'
  ENDIF 


END PROGRAM 

