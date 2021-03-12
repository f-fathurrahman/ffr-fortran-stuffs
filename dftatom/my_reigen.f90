! Solves the radial Schroedinger/Dirac eigenproblem

module my_reigen

use types, only: dp
use ode1d, only: integrate, normalize, parsefunction, get_n_nodes, get_min_idx
use my_rschroed, only: schroed_outward_adams, schroed_inward_adams
use rdirac, only: dirac_outward_adams, dirac_inward_adams
use utils, only: stop_error

implicit none
private
public solve_radial_eigenproblem, integrate_rproblem_outward

contains


subroutine integrate_rproblem_outward( l, E, R, Rp, V, Z, c, relat, &
    P, Q, imax )

integer, intent(in) :: l, relat, Z
real(dp), intent(in) :: E, c, R(:), Rp(:), V(:)
real(dp), intent(out) :: P(:), Q(:)
integer, intent(out) :: imax

integer :: kappa

if (relat == 0) then
    call schroed_outward_adams(l, Z, E, R, Rp, V, P, Q, imax)
else if (relat == 1) then
    call stop_error("Scalar relativistic case not implemented")
else if (relat == 2 .or. relat == 3) then
    if (relat == 3) then
        if (l == 0) then
              call stop_error("for l=0 only spin up (relat=2) is allowed")
        end if
        kappa = l
    else
        kappa = -l-1
    end if
    call dirac_outward_adams(c, kappa, Z, E, R, Rp, V, P, Q, imax)
else
    call stop_error("wrong value of relat.")
end if
end subroutine



!---------------------------------------------------------------
subroutine integrate_radial_problem_inward(l, E, R, Rp, V, c, &
    relat, P, Q, imin)
!---------------------------------------------------------------
  integer, intent(in) :: l, relat
  real(dp), intent(in) :: E, c, R(:), Rp(:), V(:)
  real(dp), intent(out) :: P(:), Q(:)
  integer, intent(out) :: imin
  
  integer :: kappa

  if (relat == 0) then
    call schroed_inward_adams(l, E, R, Rp, V, P, Q, imin)
  else if (relat == 1) then
    call stop_error("Scalar relativistic case not implemented")
  else if (relat == 2 .or. relat == 3) then
    if (relat == 3) then
      if (l == 0) then
        call stop_error("for l=0 only spin up (relat=2) is allowed")
      end if
      kappa = l
    else
      kappa = -l-1
    end if
    call dirac_inward_adams(c, kappa, E, R, Rp, V, P, Q, imin)
  else
    call stop_error("Wrong value of relat.")
  end if
end subroutine

!---------------------------------------------
logical function is_E_above(n, l, nods_actual)
!---------------------------------------------
    integer, intent(in) :: n, l, nods_actual
    is_E_above = nods_actual > n-l-1
end function



!-----------------------------------------------------------------------------
subroutine solve_radial_eigenproblem(n, l, Ein, eps, max_iter, &
    R, Rp, V, Z, c, relat, perturb, Emin_init, Emax_init, converged, E, P, Q)
!-----------------------------------------------------------------------------

  integer, intent(in) :: n, l, relat, Z, max_iter
  real(dp), intent(in) :: R(:), Rp(:), V(:), eps, Ein, c
  logical, intent(in) :: perturb
  real(dp), intent(in) :: Emin_init, Emax_init
  integer, intent(out) :: converged
  real(dp), intent(out) :: P(:), Q(:), E

  real(dp) :: Emin, Emax, dE, Pr(size(R)), Qr(size(R)), factor, S
  integer :: minidx, ctp, iter
  logical :: isbig
  integer :: nnodes
  logical :: last_bisect
  integer :: imin, imax

  write(*,*)
  write(*,*) 'Calling solve_radial_eigenproblem'
  write(*,*) '---------------------------------'

  E = Ein
  if(.not.(n > 0)) then
    call stop_error("n > 0 not satisfied")
  endif

  if(.not.((0 <= l) .and. (l < n))) then
    call stop_error("0 <= l < n not satisfied")
  endif

  Emax = Emax_init
  Emin = Emin_init
  if (E > Emax .or. E < Emin) E = (Emin + Emax) / 2

  iter = 0
  last_bisect = .true.
  do while (iter < max_iter)
      
    iter = iter + 1

    write(*,'(1x,A,I8,F18.10)') 'solve_reigen iter: ', iter, E
  
    ! See if bisection is converged
    if( abs(Emax - Emin) < eps ) then

      write(*,*) 'INFO: Emax and Emin difference is very small'
          
      if( .not. last_bisect ) then
        ! The perturbation theory correction was used in the last
        ! iteration and in that case, the consistent stopping criterion is
        ! to converge with abs(dE), not abs(Emax - Emin).
        ! As such we fail, because abs(Emax - Emin) is converged, but
        ! abs(dE) isn't.
        converged = 6
        return
      endif

      if( abs(Emax - Emax_init) < tiny(1._dp) ) then
        ! The algorithm didn't change Emax, so Emax_init was set incorrectly.
        write(*,*) 'Emax_init was set incorrectly'
        converged = 10
        return
      end if
          
      if( abs(Emin - Emin_init) < tiny(1._dp) ) then
        ! The algorithm didn't change Emin, so Emin_init was set incorrectly.
        write(*,*) 'Emin_init was set incorrectly'
        converged = 9
        return
      end if

      call integrate_rproblem_outward(l, E, R, Rp, V, Z, c, relat, P, Q, imax)
          
      minidx = get_min_idx(P(:imax))
      if( minidx <= 0 ) then
        ! The wavefunction doesn't have a peak
        converged = 4
        return
      end if
  
      ! Trim the wavefunction after the last minimum:
      P(minidx:) = 0.d0
      Q(minidx:) = 0.d0
  
      ! To make sure the zeros from above are not counted as nodes, we
      ! substract 1 from minidx here:
      nnodes = get_n_nodes(P(:minidx-1))
  
      if (nnodes /= n - l - 1) then
        write(*,*) 'ERROR: Wrong number of nodes for the converged energy'
        converged = 5
        return
      end if
  
      exit ! exit the while loop
    
    endif
  
    ctp = find_ctp(V + l*(l+1)/(2*R**2), E)
    
    ! If the classical turning point is too large (or cannot be found at all),
    ! we can't use inward integration to correct the energy, so we use
    ! bisection. Also use bisection if the user requests it.
    if (.not. perturb) then
      ctp = size(R)
    else if (ctp == 0) then
      ctp = size(R)
    else if (R(ctp) / R(size(R)) > 0.5_dp) then
      ctp = size(R)
    else if (size(R) - ctp <= 10) then
      ctp = size(R)
    else if (E >= 0) then
      ! Also do bisection for positive energies, as we cannot use inward
      ! integration for these
      ctp = size(R)
    end if
    
    call integrate_rproblem_outward(l, E, R(:ctp), Rp(:ctp), V(:ctp), &
      Z, c, relat, P(:ctp), Q(:ctp), imax)
    
    nnodes = get_n_nodes(P(:imax))
    write(*,*) 'nnodes = ', nnodes

    if (nnodes /= n-l-1 .or. ctp == size(R) .or. imax < ctp) then
      ! If the number of nodes is not correct, or we didn't manage to
      ! integrate all the way to "ctp", or if "ctp" was too large, we just
      ! use bisection:
      isbig = is_E_above(n, l, nnodes)
      if (isbig) then
          Emax = E
      else
          Emin = E
      end if
      E = (Emin + Emax) / 2
      last_bisect = .true.
      write(*,*) 'Need to cycle the loop'
      cycle
    end if
  
    ! Perturbation theory correction
    call integrate_radial_problem_inward( &
      l, E, R(ctp:), Rp(ctp:), V(ctp:), c, &
      relat, Pr(ctp:), Qr(ctp:), imin )

    !write(*,*) 'imin = ', imin
    if (imin > 1) then
      write(*,*) 'The inward integration didn''t integrate to the ctp'
      converged = 8
      return
    end if
  
    ! Normalize the inward solution to match the outward one:
    factor = P(ctp) / Pr(ctp)
    if (abs(factor) > 1e9) then
      write(*,*) 'Normalization factor for inward/out is too large'
      write(*,*) 'factor = ', factor
      converged = 7
      return
    end if
    Pr = Pr * factor
    Qr = Qr * factor
  
    P(ctp+1:) = Pr(ctp+1:)
    Q(ctp+1:) = Qr(ctp+1:)
    if (relat == 2 .or. relat == 3) then
      S = integrate(Rp, P**2 + Q**2)
    else
      S = integrate(Rp, P**2)
    end if
    dE = P(ctp) * (Q(ctp) - Qr(ctp)) / (2 * S)
    if (relat == 2 .or. relat == 3) then
       dE = 2 * c * dE
    end if
  
    ! The only stopping criterion for perturbation theory correction:
    if (abs(dE) < eps) exit
  
    ! We always trust the sign of dE to drive bisection
    isbig = dE < 0
    if (isbig) then
      Emax = E
    else
      Emin = E
    end if
  
    ! If the dE prediction is out of the trust region, we don't trust the value
    ! of dE, and we do bisection
    if (E + dE > Emax .or. E + dE < Emin) then
      E = (Emin + Emax) / 2
      last_bisect = .true.
    else
      E = E + dE
      last_bisect = .false.
    end if
  
  enddo

  if (iter == max_iter) then
    ! We didn't converge in 'max_iter' iterations
    converged = 2
    return
  end if
  
  ! Normalize the wavefunction:
  if (relat == 0) then
    S = integrate(Rp, P**2)
  else
    S = integrate(Rp, P**2 + Q**2)
  end if
  
  S = sqrt(abs(S))
  if (S > 0) then
    P = P / S
    Q = Q / S
  else
    ! This would happen if the function is zero, but we already check this
    ! above (converged == 4), so we fail laudly here.
    call stop_error("solve_radial_eigenproblem: zero function")
  end if
  
  converged = 0
end subroutine

!------------------------------------------
integer function find_ctp(V, E) result(ctp)
!------------------------------------------
  ! Finds the classical turning point for the potential 'V' and energy 'E'
  !
  ! Classical turning point 'ctp' is defined as E = V(ctp)
  ! The function returns the integer index into the array V.
  real(dp), intent(in) :: V(:), E
  integer :: i
  do i = size(V), 1, -1
    if (V(i)-E <= 0) then
      ctp = i
      return
    endif
  enddo
  ctp = 0
end function

end module
