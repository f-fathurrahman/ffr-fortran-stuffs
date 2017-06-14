PROGRAM test_move_alloc
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE :: a(:), b(:)

  ALLOCATE(a(3))
  a = [ 1, 2, 3 ]
  CALL move_alloc(a, b)
  PRINT *, allocated(a), allocated(b)
  PRINT *, b

END PROGRAM 

