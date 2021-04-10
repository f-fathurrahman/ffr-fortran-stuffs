subroutine my_potnucl(ptnucl, nr, r, zn, vn)

  implicit none
  ! arguments
  
  ! ptnucl : .true. if the nucleus is a point charge (in,logical)
  logical, intent(in) :: ptnucl
  
  ! nr : number of radial mesh points (in,integer)
  ! r  : radial mesh (in,real(nr))
  ! zn : nuclear charge (in,real)
  ! vn : potential on radial mesh (out,real(nr))
  integer, intent(in) :: nr
  real(8), intent(in) :: r(nr),zn
  real(8), intent(out) :: vn(nr)
  
  ! local variables
  integer ir
  real(8) rn,t1,t2
  
  ! external functions
  real(8) my_radnucl
  external my_radnucl
  
  if (zn == 0.d0) then
    vn(:) = 0.d0
    return
  endif

  if (ptnucl) then
    ! nucleus is taken to be a point particle
    vn(:) = zn/r(:)
  else
    ! approximate nuclear radius
    rn = my_radnucl(zn)
    t1 = zn/(2.d0*rn**3)
    t2 = 3.d0*rn**2
    do ir = 1,nr
      if( r(ir) < rn ) then
        vn(ir) = t1*(t2 - r(ir)**2)
      else
        vn(ir) = zn/r(ir)
      endif
    enddo
  endif
  return
end subroutine