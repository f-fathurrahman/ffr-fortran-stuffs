PROGRAM test_solve

  IMPLICIT NONE 

  ! Atomic number:
  INTEGER :: Z
  
  ! Mesh parameters:
  REAL(8), PARAMETER :: r_min = 1.0d-8, r_max = 50.d0, a = 1.0d6
  INTEGER, PARAMETER :: Nr = 3001
  REAL(8), PARAMETER :: c = 137.035999037d0, eps = 1.0d-9
  INTEGER :: n, l, relat, converged
  REAL(8), ALLOCATABLE :: R(:), u(:)
  REAL(8), ALLOCATABLE :: P(:), Q(:), Rp(:)
  REAL(8) :: Ein, E, E_exact, error
  LOGICAL, PARAMETER :: perturb = .false. !.true.

  !
  ALLOCATE( R(Nr) )
  ALLOCATE( Rp(Nr) )
  !
  ALLOCATE( u(Nr) )
  ALLOCATE( P(Nr) )
  ALLOCATE( Q(Nr) )

  CALL gen_rmesh_exp(r_min, r_max, a, Nr, R)
  CALL gen_drmesh_exp(r_min, r_max, a, Nr, Rp)
  
  ! Potential
  Z = 1
  u(:) = -Z/r

  WRITE(*,*) "Hydrogen like energies for Z = ", Z
  WRITE(*,*) "Mesh parameters (r_min, r_max, a, Nr):"
  WRITE(*,"(1x,ES10.2, F10.2, ES10.2, I10)") r_min, r_max, a, Nr
  WRITE(*,*)
  WRITE(*,"(1x,A3, A3, A15, A15, A10)") "n", "l", "E", "E_exact", "Error"
  WRITE(*,*)
  
  n = 1 ! hardcoded
  l = 0
  !do l = 0, n-1
    E_exact = - Z**2 / (2.d0 * n**2)
    Ein = -1000.d0
    relat = 0
    CALL solve_radial_eigenproblem(n, l, Ein, eps, 1000, &
        Nr, R, Rp, u, Z, c, &
        relat, perturb, -5000.0d0, 0.0d0, converged, E, P, Q)
    error = abs(E -E_exact)
    IF( converged /= 0 ) THEN
      CALL stop_error("Not converged")
    ENDIF 
    WRITE(*,"(I3, I3, F15.6, F15.6, ES10.2)") n, l, E, E_exact, error
  !end do

  WRITE(1000,*) r
  WRITE(1003,*) P
  WRITE(1004,*) Q


  DEALLOCATE( R )
  DEALLOCATE( Rp )
  DEALLOCATE( u )
  DEALLOCATE( P )
  DEALLOCATE( Q )

END PROGRAM 

