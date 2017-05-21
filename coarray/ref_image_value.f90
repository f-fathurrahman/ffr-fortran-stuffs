PROGRAM ref_image_value

  IMPLICIT NONE 

  ! XXX: does it need to be initialized?
  !INTEGER :: N[*] 
  !INTEGER, CODIMENSION[*] :: N = -99
  INTEGER :: N[*] = -99 !! cannot compile ??

  IF( this_image() == 1 ) THEN 
    WRITE(*,*) 'Program running with ', num_images(), ' images'
  ENDIF 

  N = this_image()

  SYNC ALL
  IF( this_image() == 2 ) THEN 
    WRITE(*,*) 'Before assignment:  N = ', N
    N = N[1]
    WRITE(*,*) 'After assignment:   N = ', N
  ENDIF 


END PROGRAM 

