program test_poly3
  implicit none
  real(8) :: xa(3), ya(3), x
  !
  real(8) :: poly3

  xa = (/ 1.d0, 2.d0, 3.d0 /)
  ya = xa**2 + 1.1d0

  x = 3.1d0 ! outside the function
  write(*,*) 'xa = ', xa
  write(*,*) 'ya = ', ya

  write(*,*) 'poly3 = ', poly3(xa, ya, x)
  write(*,*) 'true = ', x**2 + 1.1d0

end program