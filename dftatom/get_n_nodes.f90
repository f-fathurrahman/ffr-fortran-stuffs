SUBROUTINE get_n_nodes(N, y, nodes)
  ! Returns the number of nodes of the function 'y'
  INTEGER :: N
  REAL(8) :: y(N)
  INTEGER :: nodes
  !
  INTEGER :: last_sign, last_i, i, isy

  nodes = 0
  last_sign = int( sign(1.d0, y(1)) )
  last_i = -1

  DO i = 2,N
    isy = int( sign(1.d0, y(i)) )
    IF( isy == -last_sign ) THEN 
      last_sign = isy
      last_i = i - 1
      nodes = nodes + 1
    ENDIF 
  ENDDO 
  RETURN 
END SUBROUTINE 

