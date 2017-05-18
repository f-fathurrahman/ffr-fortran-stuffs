PROGRAM hello

  IMPLICIT NONE 
  INTEGER :: image, Nimages

  image = this_image()
  Nimages = num_images()
  
  WRITE(*,*) 'Hello, I am image: ', this_image(), ' from ', Nimages, ' Nimages'

END PROGRAM 

