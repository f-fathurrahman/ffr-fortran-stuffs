real(8) function my_radnucl(z)

  implicit none
  ! arguments
  real(8), intent(in) :: z
  
  ! local variables
  ! coefficients for computing mass number
  real(8), parameter :: c2=4.467d-3, c1=2.163d0, c0=-1.168d0
  ! coefficients for computing charge radius (fm)
  real(8), parameter :: r0=0.9071d0, r1=1.105d0, r2=-0.548d0
  ! Bohr radius in SI units (CODATA 2018)
  real(8), parameter :: br_si=0.529177210903d-10
  !
  real(8) za,a,a13,a23,a43

  za = abs(z)
  ! approximate nuclear mass number
  if (za <= 1.d0) then
    a = 1.d0
  else
    a = c2*za**2 + c1*za + c0
  endif
  
  ! approximate nuclear charge radius
  a13 = a**(1.d0/3.d0)
  a23 = a13**2
  a43 = a13*a
  my_radnucl = (r0 + r1/a23 + r2/a43)*a13
  my_radnucl = my_radnucl*1.d-15/br_si
  return
end function

