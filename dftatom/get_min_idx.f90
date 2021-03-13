SUBROUTINE get_min_idx(N, y, k)
  IMPLICIT NONE 
  ! Returns the index of the last minimum of the function 'y'
  INTEGER :: N
  REAL(8) :: y(N)
  INTEGER :: k
  
  k = N
  DO WHILE ( abs(y(k-1)) < abs(y(k) ) )
    k = k - 1
    IF( k == 1 ) EXIT 
  ENDDO 
  k = k - 1
  RETURN 
END SUBROUTINE 
