subroutine my_rschrodint(sol,l,e,nr,r,vr,nn,p0,p1,q0,q1)

  implicit none
  
  ! arguments
  
  ! sol : speed of light in atomic units (in,real)
  ! l : angular momentum quantum number (in,integer)
  ! e : energy (in,real)
  real(8), intent(in) :: sol
  integer, intent(in) :: l
  real(8), intent(in) :: e

  ! nr: number of radial mesh points (in,integer)
  ! r : radial mesh (in,real(nr))
  ! vr: potential on radial mesh (in,real(nr))
  integer, intent(in) :: nr
  real(8), intent(in) :: r(nr),vr(nr)

  ! nn: number of nodes (out,integer)
  integer, intent(out) :: nn

  ! p0: m th energy derivative of P (out,real(nr))
  ! p1: radial derivative of p0 (out,real(nr))
  ! q0: m th energy derivative of Q (out,real(nr))
  ! q1: radial derivative of q0 (out,real(nr))
  real(8), intent(out) :: p0(nr),p1(nr)
  real(8), intent(out) :: q0(nr),q1(nr)

  ! external functions
  real(8) :: poly3, poly4i

  ! local variables
  integer ir,ir0
  real(8) ri,t1,t2,t3,t4

  t1 = 1.d0/sol**2
  t2 = dble(l*(l+1))
  
  ! determine the r -> 0 boundary values of P and Q
  ri = 1.d0/r(1)
  t3 = 2.d0 + t1*(e - vr(1))
  t4 = t2/(t3*r(1)**2) + vr(1) - e
  q0(1) = 1.d0
  q1(1) = 0.d0
  p0(1) = (q1(1) + q0(1)*ri)/t4
  p1(1) = t3*q0(1) + p0(1)*ri

  ! extrapolate to the first four points
  p1(2:4) = p1(1)
  q1(2:4) = q1(1)
  nn = 0

  do ir = 2,nr
    ri = 1.d0/r(ir)
    t3 = 2.d0 + t1*(e - vr(ir))
    t4 = t2/(t3*r(ir)**2) + vr(ir) - e
    ir0 = ir-3
    if (ir0 < 1) ir0 = 1
    !
    p1(ir) = poly3( r(ir0), p1(ir0), r(ir) )
    q1(ir) = poly3( r(ir0), q1(ir0), r(ir) )
    ! integrate to find wavefunction
    p0(ir) = poly4i( r(ir0), p1(ir0), r(ir) ) + p0(ir0)
    q0(ir) = poly4i( r(ir0), q1(ir0), r(ir) ) + q0(ir0)
    ! compute the derivatives
    p1(ir) = t3*q0(ir) + p0(ir)*ri
    q1(ir) = t4*p0(ir) - q0(ir)*ri
    ! integrate for correction
    p0(ir) = poly4i( r(ir0), p1(ir0), r(ir) ) + p0(ir0)
    q0(ir) = poly4i( r(ir0), q1(ir0), r(ir) ) + q0(ir0)
    ! compute the derivatives again
    p1(ir) = t3*q0(ir) + p0(ir)*ri
    q1(ir) = t4*p0(ir) - q0(ir)*ri
    ! check for overflow
    if ((abs(p0(ir)) > 1.d100) .or. ( abs(p1(ir)) > 1.d100) .or. &
        (abs(q0(ir)) > 1.d100) .or. ( abs(q1(ir)) > 1.d100)) then
      p0(ir:nr) = p0(ir)
      p1(ir:nr) = p1(ir)
      q0(ir:nr) = q0(ir)
      q1(ir:nr) = q1(ir)
      return
    end if
    ! check for node
    if( p0(ir-1)*p0(ir) < 0.d0) nn = nn + 1
  end do
  return
end subroutine
