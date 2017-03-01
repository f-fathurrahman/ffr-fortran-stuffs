PROGRAM linked_list
  IMPLICIT NONE 
  REAL(8), TARGET :: stat1
  REAL(8), TARGET :: stat2
  REAL(8), POINTER :: ptr1

  stat1 = 1.1d0
  stat2 = 2.3d0

  ptr1 => stat1
  WRITE(*,*) 'ptr1 = ', ptr1

  ptr1 => stat2
  WRITE(*,*) 'ptr1 = ', ptr1
END PROGRAM 

