!-----------------------------------------------------------------------------
SUBROUTINE solve_radial_eigenproblem(n, l, Ein, eps, max_iter, &
    Nr, R, Rp, V, Z, c, relat, perturb, Emin_init, Emax_init, converged, E, P, Q)
!-----------------------------------------------------------------------------

  IMPLICIT NONE 
  INTEGER :: Nr
  INTEGER, INTENT(in) :: n, l, relat, Z, max_iter
  REAL(8), INTENT(in) :: R(Nr), Rp(Nr), V(Nr), eps, Ein, c
  LOGICAL, INTENT(in) :: perturb
  REAL(8), INTENT(in) :: Emin_init, Emax_init
  INTEGER, INTENT(out) :: converged
  REAL(8), INTENT(out) :: P(Nr), Q(Nr), E

  REAL(8) :: Emin, Emax, dE, factor, S
  REAL(8), ALLOCATABLE :: Pr(:), Qr(:)
  INTEGER :: minidx, ctp, iter
  LOGICAL :: isbig
  INTEGER :: nnodes
  LOGICAL :: last_bisect
  INTEGER :: imin, imax

  WRITE(*,*)
  WRITE(*,*) 'Calling solve_radial_eigenproblem'
  WRITE(*,*) '---------------------------------'

  ALLOCATE( Pr(Nr) )
  ALLOCATE( Qr(Nr) )

  E = Ein
  IF( .not. (n > 0) ) THEN
    CALL stop_error("n > 0 not satisfied")
  ENDIF

  IF( .not. ( (0 <= l) .and. (l < n) ) ) THEN 
    CALL stop_error("0 <= l < n not satisfied")
  ENDIF

  Emax = Emax_init
  Emin = Emin_init
  IF( E > Emax .or. E < Emin ) THEN 
    E = (Emin + Emax) / 2
  ENDIF 

  iter = 0
  last_bisect = .true.

  DO WHILE( iter < max_iter )
      
    iter = iter + 1

    WRITE(*,*)
    WRITE(*,'(1x,A,I8,F18.10)') 'solve_reigen iter: ', iter, E
  
    ! See if bisection is converged
    IF( abs(Emax - Emin) < eps ) THEN

      WRITE(*,*) 'INFO: Emax and Emin difference is very small'
          
      IF( .not. last_bisect ) THEN 
        ! The perturbation theory correction was used in the last
        ! iteration and in that case, the consistent stopping criterion is
        ! to converge with abs(dE), not abs(Emax - Emin).
        ! As such we fail, because abs(Emax - Emin) is converged, but
        ! abs(dE) isn't.
        WRITE(*,*) 'ERROR: abs(dE) is not converged'
        converged = 6
        RETURN 
      ENDIF 

      IF( abs(Emax - Emax_init) < tiny(1.d0) ) THEN 
        WRITE(*,*) "ERROR: The algorithm didn't change Emax"
        WRITE(*,*) 'ERROR: Emax_init was set incorrectly'
        converged = 10
        RETURN
      ENDIF 
          
      IF( abs(Emin - Emin_init) < tiny(1.d0) ) THEN 
        WRITE(*,*) "The algorithm didn't change Emin"
        WRITE(*,*) 'Emin_init was set incorrectly'
        converged = 9
        RETURN 
      ENDIF 

      WRITE(*,*) 'DEBUG: integrate_rproblem_outward until Nr'
      CALL integrate_rproblem_outward( Nr, l, E, R, Rp, V, Z, c, relat, P, Q, imax )
      
      CALL get_min_idx( imax, P(:imax), minidx )
      WRITE(*,*) 'minidx = ', minidx

      IF( minidx <= 0 ) THEN 
        WRITE(*,*) "ERROR: The wavefunction doesn't have a peak"
        converged = 4
        RETURN 
      ENDIF 
  
      ! Trim the wavefunction after the last minimum:
      P(minidx:) = 0.d0
      Q(minidx:) = 0.d0
  
      ! To make sure the zeros from above are not counted as nodes, we
      ! substract 1 from minidx here:
      CALL get_n_nodes( minidx-1, P(:minidx-1), nnodes )
      WRITE(*,*) 'DEBUG: nnodes = ', nnodes
  
      IF( nnodes /= n - l - 1) THEN 
        WRITE(*,*) 'ERROR: Wrong number of nodes for the converged energy'
        converged = 5
        RETURN 
      ENDIF 

      WRITE(*,*) 'DEBUG: will exit the while loop'
      EXIT ! exit the while loop
    
    ENDIF 
  
    CALL find_ctp( Nr, V + l*(l+1)/(2*R**2), E, ctp )
    WRITE(*,*) 'DEBUG: find_ctp, ctp = ', ctp

    ! If the classical turning point is too large (or cannot be found at all),
    ! we can't use inward integration to correct the energy, so we use
    ! bisection. Also use bisection if the user requests it.
    IF( .not. perturb ) THEN 
      WRITE(*,*) 'DEBUG: not using perturbation correction, ctp set to Nr'
      ctp = Nr
    !
    ELSEIF( ctp == 0 ) THEN 
      ctp = Nr
    !
    ELSEIF( R(ctp) / R(Nr) > 0.5d0) THEN 
      ctp = Nr
    !
    ELSEIF( Nr - ctp <= 10 ) THEN 
      ctp = Nr
    !
    ELSEIF( E >= 0 ) THEN 
      ! Also do bisection for positive energies, as we cannot use inward
      ! integration for these
      ctp = Nr
    ENDIF 
    
    WRITE(*,*) 'DEBUG: integrate the problem outward until ctp = ', ctp
    CALL integrate_rproblem_outward( ctp, l, E, R(:ctp), Rp(:ctp), V(:ctp), &
      Z, c, relat, P(:ctp), Q(:ctp), imax )

    WRITE(*,*) 'imax = ', imax

    CALL get_n_nodes( imax, P(:imax), nnodes )
    WRITE(*,*) 'DEBUG: nnodes = ', nnodes


    ! If the number of nodes is not correct, or we didn't manage to
    ! integrate all the way to "ctp", or if "ctp" was too large, we just
    ! use bisection:
    IF( nnodes /= n-l-1 .or. ctp == Nr .or. imax < ctp ) THEN 

      WRITE(*,'(1x,A,2F18.10)') 'DEBUG: before, Emin, Emax = ', Emin, Emax

      CALL is_E_above(n, l, nnodes, isbig)
      !
      WRITE(*,*) 'DEBUG: isbig = ', isbig
      !
      IF( isbig ) THEN 
        Emax = E
      ELSE 
        Emin = E
      ENDIF 
      
      WRITE(*,'(1x,A,2F18.10)') 'DEBUG: after , Emin, Emax = ', Emin, Emax
      
      WRITE(*,*) 'DEBUG: Updating E via bisection'
      E = (Emin + Emax) / 2
      last_bisect = .true.
      WRITE(*,*) 'Need to cycle the loop'
      !STOP 'ffr'
      CYCLE
    ENDIF 
  
    ! Perturbation theory correction
    WRITE(*,*) 'DEBUG: Perturbation theory correction'
    WRITE(*,*) 'DEBUG: perturb = ', perturb
    CALL integrate_rproblem_inward( &
      Nr-ctp+1, l, E, R(ctp:), Rp(ctp:), V(ctp:), c, &
      relat, Pr(ctp:), Qr(ctp:), imin )

    !write(*,*) 'imin = ', imin
    IF( imin > 1 ) THEN 
      WRITE(*,*) "ERROR: The inward integration didn't integrate to the ctp"
      converged = 8
      RETURN 
    ENDIF 
  
    ! Normalize the inward solution to match the outward one:
    factor = P(ctp) / Pr(ctp)
    IF( abs(factor) > 1e9 ) THEN 
      WRITE(*,*) 'ERROR: Normalization factor for inward/out is too large'
      WRITE(*,*) 'ERROR: factor = ', factor
      converged = 7
      RETURN 
    ENDIF 
    Pr = Pr * factor
    Qr = Qr * factor
  
    P(ctp+1:) = Pr(ctp+1:)
    Q(ctp+1:) = Qr(ctp+1:)
    IF( relat == 2 .or. relat == 3 ) THEN 
      CALL integrate_trapz_7( Nr, Rp, P**2 + Q**2, S )
    ELSE 
      CALL integrate_trapz_7( Nr, Rp, P**2, S )
    ENDIF 
    dE = P(ctp) * (Q(ctp) - Qr(ctp)) / (2 * S)
    
    IF( relat == 2 .or. relat == 3 ) THEN 
       dE = 2*c*dE
    ENDIF 
  
    ! The only stopping criterion for perturbation theory correction:
    if (abs(dE) < eps) EXIT
  
    ! We always trust the sign of dE to drive bisection
    isbig = dE < 0
    IF( isbig ) THEN 
      Emax = E
    ELSE 
      Emin = E
    ENDIF 
  
    ! If the dE prediction is out of the trust region, we don't trust the value
    ! of dE, and we do bisection
    IF( E + dE > Emax .or. E + dE < Emin ) THEN 
      E = (Emin + Emax) / 2
      last_bisect = .true.
    ELSE 
      E = E + dE
      last_bisect = .false.
    ENDIF 
  
  ENDDO ! end of bisection iterations

  IF( iter == max_iter ) THEN 
    WRITE(*,*) "WARNING: We didn't converge after max_iter", max_iter
    converged = 2
    RETURN 
  ENDIF 
  
  ! Normalize the wavefunction:
  IF( relat == 0 ) THEN 
    CALL integrate_trapz_7(Nr, Rp, P**2, S)
  ELSE 
    CALL integrate_trapz_7(Nr, Rp, P**2 + Q**2, S)
  ENDIF 
  
  S = sqrt(abs(S))
  IF( S > 0 ) THEN 
    P = P / S
    Q = Q / S
  ELSE 
    ! This would happen if the function is zero, but we already check this
    ! above (converged == 4), so we fail laudly here.
    CALL stop_error("solve_radial_eigenproblem: zero function")
  ENDIF
  
  converged = 0

  DEALLOCATE( Pr )
  DEALLOCATE( Qr )

  RETURN 
END SUBROUTINE 

