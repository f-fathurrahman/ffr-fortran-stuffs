program test_poly4i
  implicit none
  real(8) :: xa(4), ya(4), x
  !
  real(8) :: poly4i

  xa = (/ 1.1d0, 2.d0, 3.d0, 4.d0 /)
  ya = xa**2 + 1.1d0

  x = 4.1d0
  write(*,*) 'xa = ', xa
  write(*,*) 'ya = ', ya

  write(*,*) 'poly4i = ', poly4i(xa, ya, x)
  !write(*,*) 'true = ', x**2 + 1.1d0

end program