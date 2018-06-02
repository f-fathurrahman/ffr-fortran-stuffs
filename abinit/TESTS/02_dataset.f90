PROGRAM test
  USE defs_abitypes
  IMPLICIT NONE 
  TYPE(dataset_type) :: dtset

  dtset%usewvl = 0

  WRITE(*,*) dtset%iomode
  WRITE(*,*) dtset%usewvl
  WRITE(*,*) dtset%prtgsr
  
  WRITE(*,*)
  WRITE(*,*) 'Pass here ...'
END PROGRAM 

