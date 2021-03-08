subroutine my_rdiracint(sol, kpa, e, nr, r, vr, nn, g0, g1, f0, f1)

  implicit none
  
  ! arguments

  ! sol  : speed of light in atomic units (in,real)
  ! kpa  : quantum number kappa (in,integer)
  ! e    : energy (in,real)
  ! nr   : number of radial mesh points (in,integer)
  ! r    : radial mesh (in,real(nr))
  ! vr   : potential on radial mesh (in,real(nr))
  real(8), intent(in) :: sol
  integer, intent(in) :: kpa
  real(8), intent(in) :: e
  integer, intent(in) :: nr
  real(8), intent(in) :: r(nr), vr(nr)

  ! nn   : number of nodes (out,integer)
  integer, intent(out) :: nn

  ! g0   : m th energy derivative of the major component multiplied by r
  !        (out,real(nr))
  ! g1   : radial derivative of g0 (out,real(nr))
  ! f0   : m th energy derivative of the minor component multiplied by r
  !        (out,real(nr))
  ! f1   : radial derivative of f0 (out,real(nr))
  real(8), intent(out) :: g0(nr),g1(nr)
  real(8), intent(out) :: f0(nr),f1(nr)

  ! external functions
  real(8) :: poly3, poly4i

  ! local variables
  integer ir,ir0
  ! rescaling limit
  real(8), parameter :: rsc=1.d100
  real(8) ci,e0,t1,t2,t3,t4

  if (nr < 4) then
    write(*,*)
    write(*,'("Error(my_rdiracint): nr < 4 : ",I8)') nr
    write(*,*)
    stop
  end if

  ! inverse speed of light
  ci = 1.d0/sol
  
  ! electron rest energy
  e0 = sol**2
  t1 = 2.d0*e0 + e
  
  ! determine the r -> 0 boundary values of F and G
  t2 = dble(kpa)/r(1)
  t3 = ci*(t1 - vr(1))
  t4 = ci*(vr(1) - e)
  f0(1) = 1.d0
  f1(1) = 0.d0
  g0(1) = (f1(1) - t2*f0(1))/t4
  g1(1) = t3*f0(1) - t2*g0(1)

  ! extrapolate to the first four points
  g1(2:4) = g1(1)
  f1(2:4) = f1(1)
  
  nn = 0
  do ir = 2,nr
    t2 = dble(kpa)/r(ir)
    t3 = ci*(t1 - vr(ir))
    t4 = ci*(vr(ir) - e)
    ir0 = ir - 3
    if (ir0 < 1) ir0 = 1
    !
    g1(ir) = poly3( r(ir0), g1(ir0), r(ir) )
    f1(ir) = poly3( r(ir0), f1(ir0), r(ir) )
    ! integrate to find wavefunction
    g0(ir) = poly4i( r(ir0), g1(ir0), r(ir) ) + g0(ir0)
    f0(ir) = poly4i( r(ir0), f1(ir0), r(ir) ) + f0(ir0)
    ! compute the derivatives
    g1(ir) = t3*f0(ir) - t2*g0(ir)
    f1(ir) = t4*g0(ir) + t2*f0(ir)
    ! integrate for correction
    g0(ir) = poly4i( r(ir0), g1(ir0), r(ir) ) + g0(ir0)
    f0(ir) = poly4i( r(ir0), f1(ir0), r(ir) ) + f0(ir0)
    ! compute the derivatives again
    g1(ir) = t3*f0(ir) - t2*g0(ir)
    f1(ir) = t4*g0(ir) + t2*f0(ir)
    ! check for overflow
    if ((abs(g0(ir)) > rsc) .or. (abs(g1(ir)) > rsc) .or. &
        (abs(f0(ir)) > rsc) .or. (abs(f1(ir)) > rsc)) then
      ! set the remaining points and return
      g0(ir:nr) = g0(ir)
      g1(ir:nr) = g1(ir)
      f0(ir:nr) = f0(ir)
      f1(ir:nr) = f1(ir)
      return
    end if
    ! check for node
    if( g0(ir-1)*g0(ir) < 0.d0 ) nn = nn + 1
  end do
  return

end subroutine

