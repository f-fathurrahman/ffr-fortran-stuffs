PROGRAM cobound

  IMPLICIT NONE 
  CHARACTER(10) :: i[-2:2,2,1:*]

  IF( this_image() == num_images() ) THEN 
    WRITE(*,*) 'this_image = ', this_image()
    WRITE(*,*) 'this_image(i) = ', this_image(i)
    WRITE(*,*) 'lcobound(i) = ', lcobound(i)
    WRITE(*,*) 'ucobound(i) = ', ucobound(i)
  ENDIF 

END PROGRAM 

