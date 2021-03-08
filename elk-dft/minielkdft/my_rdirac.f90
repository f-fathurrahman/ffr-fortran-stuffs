subroutine my_rdirac( sol, n, l, k, nr, r, vr, eval, g0, f0)

  implicit none
  ! arguments
  
  ! sol  : speed of light in atomic units (in,real)
  real(8), intent(in) :: sol

  ! n  : principal quantum number (in,integer)
  ! l  : quantum number l (in,integer)
  ! k  : quantum number k (l or l+1) (in,integer)
  ! nr : number of radial mesh points (in,integer)
  integer, intent(in) :: n,l,k,nr

  ! r    : radial mesh (in,real(nr))
  ! vr   : potential on radial mesh (in,real(nr))
  real(8), intent(in) :: r(nr),vr(nr)

  ! eval : eigenvalue without rest-mass energy (inout,real)
  real(8), intent(inout) :: eval

  ! g0   : major component of the radial wavefunction (out,real(nr))
  ! f0   : minor component of the radial wavefunction (out,real(nr))
  real(8), intent(out) :: g0(nr),f0(nr)
  
  ! local variables
  integer, parameter :: maxit=2000

  integer kpa,it,ir,irm
  integer nn,nnd,nndp
  
  ! energy convergence tolerance
  real(8), parameter :: eps=1.d-12
  real(8) t1,de
  
  ! automatic arrays
  real(8) g1(nr),f1(nr),fr(nr)
  
  ! external functions
  real(8) splint
  external splint


  if (k <= 0) then
    write(*,*)
    write(*,'("Error(my_rdirac): k <= 0 : ",I8)') k
    write(*,*)
    stop
  end if

  if (k > n) then
    write(*,*)
    write(*,'("Error(my_rdirac): incompatible n and k : ",2I8)') n,k
    write(*,*)
    stop
  end if

  if ((k == n) .and. (l /= k-1)) then
    write(*,*)
    write(*,'("Error(my_rdirac): incompatible n, k and l : ",3I8)') n,k,l
    write(*,*)
    stop
  end if

  if (k == l) then
    kpa = k
  else if (k == l+1) then
    kpa = -k
  else
    write(*,*)
    write(*,'("Error(my_rdirac): incompatible l and k : ",2I8)') l,k
    write(*,*)
    stop
  endif

  if (nr < 4) then
    write(*,*)
    write(*,'("Error(my_rdirac): nr < 4 : ",I8)') nr
    write(*,*)
    stop
  end if

  de = 1.d0
  nndp = 0

  do it = 1,maxit
    
    ! integrate the Dirac equation
    call my_rdiracint(sol, kpa, eval, nr, r, vr, nn, g0, g1, f0, f1)
    
    ! check the number of nodes
    nnd = nn - (n-l-1)
    if (nnd > 0) then
      eval = eval - de
    else
      eval = eval + de
    end if
    
    if (it > 1) then
      if ((nnd /= 0) .or. (nndp /= 0)) then
        if (nnd*nndp <= 0) then
          de = de*0.5d0
        else
          de = de*1.1d0
        endif
      endif
    end if
    
    nndp = nnd
    if( de < eps*(abs(eval) + 1.d0) ) goto 20
  end do
  
  write(*,*)
  write(*,'("Warning(rdirac): maximum iterations exceeded")')
  20 continue
  
  ! find effective infinity and set wavefunction to zero after that point
  
  ! major component
  irm = nr
  do ir = 2,nr
    if( (g0(ir-1)*g0(ir) < 0.d0) .or. (g1(ir-1)*g1(ir) < 0.d0) ) irm = ir
  enddo
  g0(irm:nr) = 0.d0
  
  ! minor component
  irm = nr
  do ir = 2,nr
    if( (f0(ir-1)*f0(ir) < 0.d0) .or. (f1(ir-1)*f1(ir) < 0.d0) ) irm = ir
  end do
  f0(irm:nr) = 0.d0
  
  ! normalise
  do ir = 1,nr
    fr(ir) = g0(ir)**2 + f0(ir)**2
  end do
  t1 = splint(nr,r,fr)
  t1 = sqrt(abs(t1))
  !
  if (t1 > 0.d0) then
    t1 = 1.d0/t1
  else
    write(*,*)
    write(*,'("Error(my_rdirac): zero wavefunction")')
    write(*,*)
    stop
  end if
  !
  call dscal(nr, t1, g0, 1)
  call dscal(nr, t1, f0, 1)
  return
end subroutine
