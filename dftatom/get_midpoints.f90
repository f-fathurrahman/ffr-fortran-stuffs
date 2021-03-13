SUBROUTINE get_midpoints(N, R, V, Vmid)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8), INTENT(in) :: R(N), V(N)
  REAL(8) :: Vmid(N-1)
  INTEGER :: i

  IF( .not. (size(R) == size(V)) ) THEN 
      CALL stop_error("get_midpoints: incorrect array sizes")
  ENDIF 

  DO i = 1, N - 1
    call get_interp_val(V, (R(i) + R(i+1))/2, R, N, Vmid(i), i+1)
  ENDDO 
END SUBROUTINE  

