MODULE module1
  INTEGER :: var1, var2, N
  REAL(8), ALLOCATABLE :: arr1(:)
END

SUBROUTINE init_arr()
  USE module1, ONLY : N, arr1, var1, var2
  IMPLICIT NONE
  INTEGER :: i

  N    = 5  ! set 
  var1 = N*3 - 1
  var2 = var1 - 2
  !
  ALLOCATE( arr1(N) )
  DO i = 1, N
    CALL random_number( arr1(i) )
  ENDDO
END SUBROUTINE

SUBROUTINE fin_arr()
  USE module1, ONLY : arr1
  DEALLOCATE( arr1 )
END SUBROUTINE

PROGRAM t_module
  USE module1, ONLY : Nelem=>N, arr1, var1, var2
  IMPLICIT NONE
  INTEGER :: i

  WRITE(*,'(A,3I4)') 'var1, var2, Nelem = ', var1, var2, Nelem
  CALL init_arr()
  WRITE(*,'(A,3I4)') 'var1, var2, Nelem = ', var1, var2, Nelem

  WRITE(*,*) 'arr1:'
  DO i = 1, Nelem
    WRITE(*,'(I5,F18.10)') i, arr1(i)
  ENDDO
  
  CALL fin_arr()
END PROGRAM
