!---------------------------------------------------
REAL(8) FUNCTION integrate_trapz_1(Rp, f) RESULT(s)
!---------------------------------------------------
  IMPLICIT NONE
  REAL(8), INTENT(in) :: Rp(:), f(:)
  
  REAL(8) :: g(size(Rp))
  INTEGER :: N
  N = size(Rp)
  g = f * Rp
  s = (g(1) + g(N)) / 2
  s = s + sum(g(2:N-1))
END FUNCTION 

!---------------------------------------------------
REAL(8) FUNCTION integrate_trapz_3(Rp, f) RESULT(s)
!---------------------------------------------------
  IMPLICIT NONE
  REAL(8), INTENT(in) :: Rp(:), f(:)
  
  REAL(8) :: g(size(Rp))
  INTEGER :: N
  N = size(Rp)
  g = f * Rp
  s = (9 * (g(1) + g(N)) + 28 * (g(2) + g(N-1)) + 23 * (g(3) + g(N-2))) / 24
  s = s + sum(g(4:N-3))
END FUNCTION 

!---------------------------------------------------
REAL(8) FUNCTION integrate_trapz_5(Rp, f) RESULT(s)
!---------------------------------------------------
  IMPLICIT NONE
  REAL(8), INTENT(in) :: Rp(:), f(:)
  
  REAL(8) :: g(size(Rp))
  INTEGER :: N
  N = size(Rp)
  g = f * Rp
  s = (  475 * (g(1) + g(N  )) &
      + 1902 * (g(2) + g(N-1)) &
      + 1104 * (g(3) + g(N-2)) &
      + 1586 * (g(4) + g(N-3)) &
      + 1413 * (g(5) + g(N-4)) &
      ) / 1440
  s = s + sum(g(6:N-5))
END FUNCTION 

!---------------------------------------------------
SUBROUTINE integrate_trapz_7(N, Rp, f, s)
!---------------------------------------------------
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8), INTENT(in) :: Rp(N), f(N)
  REAL(8), INTENT(out) :: s
  !
  REAL(8), ALLOCATABLE :: g(:) ! automatic array ?

  ALLOCATE( g(N) )

  g = f * Rp
  s = (  36799 * (g(1) + g(N  )) &
      + 176648 * (g(2) + g(N-1)) &
      +  54851 * (g(3) + g(N-2)) &
      + 177984 * (g(4) + g(N-3)) &
      +  89437 * (g(5) + g(N-4)) &
      + 130936 * (g(6) + g(N-5)) &
      + 119585 * (g(7) + g(N-6)) &
      ) / 120960
  s = s + sum(g(8:N-7))

  DEALLOCATE( g )
  RETURN 
END SUBROUTINE 

