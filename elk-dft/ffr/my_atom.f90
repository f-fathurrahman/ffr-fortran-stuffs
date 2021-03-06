subroutine my_atom(sol, ptnucl, &
    & zn, nst, n, l, k, occ, &
    & xctype, xcgrad, &
    & nr, r, &
    & eval, rho, vr, rwf )

  use modxcifc, only: xcifc

  implicit none
  ! arguments
  
  ! speed of light in atomic units (in,real)
  real(8), intent(in) :: sol

  ! .true. if the nucleus is a point particle (in,logical)
  logical, intent(in) :: ptnucl

  ! nuclear charge (in,real)
  real(8), intent(in) :: zn
  
  ! number of states to solve for (in,integer)
  integer, intent(in) :: nst

  ! n   : priciple quantum number of each state (in,integer(nst))
  ! l   : quantum number l of each state (in,integer(nst))
  ! k   : quantum number k (l or l+1) of each state (in,integer(nst))
  ! occ : occupancy of each state (inout,real(nst)) ! FIXME: why inout?
  integer, intent(in) :: n(nst), l(nst), k(nst)
  real(8), intent(inout) :: occ(nst)

  ! xctype : exchange-correlation type (in,integer(3))
  ! xcgrad : 1 for GGA functional, 0 otherwise (in,integer)
  integer, intent(in) :: xctype(3), xcgrad

  ! nr : number of radial mesh points (in,integer)
  ! r  : radial mesh (in,real(nr))
  integer, intent(in) :: nr
  real(8), intent(in) :: r(nr)

  !   eval   : eigenvalue without rest-mass energy for each state (out,real(nst))
  !   rho    : charge density (out,real(nr))
  !   vr     : self-constistent potential (out,real(nr))
  !   rwf    : major and minor components of radial wavefunctions for each state
  !            (out,real(nr,2,nst))
  real(8), intent(out) :: eval(nst)
  real(8), intent(out) :: rho(nr), vr(nr)
  real(8), intent(out) :: rwf(nr,2,nst)
  
  ! local variables
  integer, parameter :: maxscl=200
  integer ir,ist,iscl
  real(8), parameter :: fourpi=12.566370614359172954d0
  
  ! potential convergence tolerance
  real(8), parameter :: eps=1.d-6
  real(8) ss, dv, dvp, ze, beta, t1
  
  ! allocatable arrays
  real(8), allocatable :: vn(:), vh(:), ex(:), ec(:), vx(:), vc(:), vrp(:)
  real(8), allocatable :: ri(:), wpr(:,:), fr1(:), fr2(:), gr1(:), gr2(:)
  real(8), allocatable :: grho(:), g2rho(:), g3rho(:)

  if (nst.le.0) then
    write(*,*)
    write(*,'("Error(atom): invalid nst : ",I8)') nst
    write(*,*)
    stop
  end if

  write(*,*) '--- Solving Dirac-Kohn-Sham equation for atom ---'
  write(*,*) 'zn     = ', zn
  write(*,*) 'nst    = ', nst
  write(*,*) 'sol    = ', sol
  write(*,*) 'ptnucl = ', ptnucl

  do ist = 1,nst
    write(*,'(1x,A,I5)', advance="no") ' n = ', n(ist)
    write(*,'(1x,A,I5)', advance="no") ' l = ', l(ist)
    write(*,'(1x,A,I5)', advance="no") ' k = ', k(ist)
    write(*,*)
 enddo

  ! allocate local arrays
  allocate(vn(nr), vh(nr), ex(nr), ec(nr), vx(nr), vc(nr), vrp(nr))
  allocate(ri(nr), wpr(4,nr), fr1(nr), fr2(nr), gr1(nr), gr2(nr))
  if (xcgrad == 1) then
    allocate(grho(nr), g2rho(nr), g3rho(nr))
  end if

  ! find total electronic charge
  ze = 0.d0
  do ist = 1,nst
    ze = ze + occ(ist)
  enddo

  ! set up nuclear potential
  call potnucl(ptnucl, nr, r, zn, vn)
  do ir = 1,nr
    ri(ir) = 1.d0/r(ir)
    ! initialise the Kohn-Sham potential to the nuclear potential
    vr(ir) = vn(ir)
  enddo

  ! determine the weights for radial integration
  call wsplintp(nr, r, wpr)
  
  dvp = 0.d0
  vrp(:) = 0.d0
  
  ! initialise mixing parameter
  beta = 0.5d0
  
  ! initialise eigenvalues to relativistic values (minus the rest mass energy)
  do ist = 1,nst
    t1 = sqrt(dble(k(ist)**2) - (zn/sol)**2)
    t1 = (dble(n(ist) - abs(k(ist))) + t1)**2
    t1 = 1.d0 + ((zn/sol)**2)/t1
    eval(ist) = sol**2/sqrt(t1) - sol**2
  enddo
  
  ! start of self-consistent loop
  do iscl = 1,maxscl
    
    ! solve the Dirac equation for each state
    do ist = 1,nst
      call rdirac( sol, n(ist), l(ist), k(ist), &
                 & nr, r, vr, &
                 & eval(ist), rwf(:,1,ist), rwf(:,2,ist) )
    enddo

    ! compute the charge density
    do ir = 1,nr
      sum = 0.d0
      do ist = 1,nst
        sum = sum + occ(ist)*(rwf(ir,1,ist)**2 + rwf(ir,2,ist)**2)
      enddo
      fr1(ir) = sum
      fr2(ir) = sum*ri(ir)
      rho(ir) = (1.d0/fourpi)*sum*ri(ir)**2
    enddo
    call splintwp(nr, wpr, fr1, gr1)
    call splintwp(nr, wpr, fr2, gr2)
    
    ! find the Hartree potential
    t1 = gr2(nr)
    do ir = 1,nr
      vh(ir) = gr1(ir)*ri(ir) + t1 - gr2(ir)
    end do
    
    ! normalise charge density and potential
    t1 = ze/gr1(nr)
    rho(:) = t1*rho(:)
    vh(:) = t1*vh(:)
    
    ! compute the exchange-correlation energy and potential
    if (xcgrad == 1) then
      ! GGA functional
      ! |grad rho|
      call fderiv(1, nr, r, rho, grho)
      ! grad^2 rho
      call fderiv(2, nr, r, rho, g2rho)
      do ir=1,nr
        g2rho(ir) = g2rho(ir) + 2.d0*ri(ir)*grho(ir)
      enddo
      ! approximate (grad rho).(grad |grad rho|)
      do ir = 1,nr
        g3rho(ir) = grho(ir)*g2rho(ir)
      end do
      call xcifc(xctype, n=nr, rho=rho, grho=grho, g2rho=g2rho, g3rho=g3rho, &
              &  ex=ex, ec=ec, vx=vx, vc=vc)
    else
      ! LDA functional
      call xcifc(xctype, n=nr, rho=rho, ex=ex, ec=ec, vx=vx, vc=vc)
    endif
    
    ! self-consistent potential
    vr(:) = vh(:) + vx(:) + vc(:)
    
    ! determine change in potential
    ss = 0.d0
    do ir = 1,nr
      ss = ss + ( vr(ir) - vrp(ir) )**2
    enddo
    dv = sqrt(ss)/dble(nr)
    
    if (iscl > 2) then
      ! reduce beta if change in potential is diverging
      if (dv > dvp) beta = beta*0.8d0
      beta = max(beta, 0.01d0)
    endif
    dvp = dv

    do ir = 1,nr
      ! mix old and new potentials
      vr(ir) = (1.d0 - beta)*vrp(ir) + beta*vr(ir)
      vrp(ir) = vr(ir) ! old potential ? ffr
      ! add nuclear potential
      vr(ir) = vr(ir) + vn(ir)
    enddo
    
    ! check for convergence
    if ( (iscl > 2) .and. (dv < eps) ) goto 10
    
  ! end self-consistent loop
  enddo
  
  write(*,*)
  write(*,'("Warning(atom): maximum iterations exceeded")')
  
  10 continue
  deallocate(vn,vh,ex,ec,vx,vc,vrp)
  deallocate(ri,wpr,fr1,fr2,gr1,gr2)
  if (xcgrad == 1) deallocate(grho,g2rho,g3rho)
  return
end subroutine

