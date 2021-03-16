SUBROUTINE V_to_rho( d )
  USE m_dft_data, ONLY: dft_data_t
  IMPLICIT NONE
  ! Calculates rho from V by solving Kohn-Sham equations
  TYPE(dft_data_t), INTENT(inout) :: d

  REAL(8), ALLOCATABLE :: P(:), Q(:), Y(:)
  INTEGER :: converged, i, n, l, relat
  REAL(8) :: Ein, Emin_init, Emax_init
  INTEGER :: Nr, n_orb
  REAL(8), PARAMETER :: pi=4.d0*atan(1.d0)

  Nr = size(d%R)
  n_orb = size(d%no)

  ALLOCATE( P(Nr) )
  ALLOCATE( Q(Nr) )
  ALLOCATE( Y(Nr) )

  d%rho(:) = 0.d0

  !print *, d%scf_iter, d%ks_energies
  DO i = 1,n_orb
    !
    n = d%no(i)
    l = d%lo(i)
    !
    IF( d%dirac ) THEN 
      IF( d%so(i) == 1) then
        relat = 2
      ELSE 
        relat = 3
      ENDIF 
    ELSE 
      relat = 0
    ENDIF 
  
    Ein = d%ks_energies(i)
    Emax_init = d%Emax_init(i)
    Emin_init = d%Emin_init(i)

    CALL solve_radial_eigenproblem( n, l, Ein, d%reigen_eps, &
        d%reigen_max_iter, &
        Nr, d%R, d%Rp, d%V_tot, &
        d%Z, d%c, relat, d%perturb, Emin_init, Emax_init, &
        converged, d%ks_energies(i), P, Q )
    
    IF( converged /= 0 ) THEN 
      WRITE(*,*) "converged=", converged
      WRITE(*,*) d%scf_iter, n, l, relat
      !print *, "skipping the state"
      !Y = 0
      CALL stop_error("V_to_rho: Radial eigen problem didn't converge")
    ENDIF 

    IF( relat == 0 ) THEN 
      Y = P / d%R
    ELSE 
      Y = sqrt(P**2 + Q**2) / d%R
    ENDIF 
    d%rho = d%rho + d%fo(i) * Y**2
    d%orbitals(:, i) = Y
  ENDDO 
  d%rho = d%rho / (4*pi)

  RETURN 
END SUBROUTINE 

