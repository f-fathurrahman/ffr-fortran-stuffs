PROGRAM test

  IMPLICIT NONE 
  REAL(8) :: r_min, r_max, a
  INTEGER :: Nr
  REAL(8), ALLOCATABLE :: rmesh(:), drmesh(:)
  INTEGER :: i
  REAL(8), ALLOCATABLE :: V(:), Vmid(:)
  REAL(8), ALLOCATABLE :: u1(:), u2(:)
  INTEGER :: Z, l, n
  INTEGER :: iterShoot
  INTEGER :: imax, ctp
  REAL(8) :: E, Emin, Emax, dE
  REAL(8) :: S
  INTEGER :: minidx, nnodes
  LOGICAL :: isbig

  ! Setup radial mesh
  Nr = 3001
  r_min = 1.d-9
  r_max = 50.d0
  a = 1d6
  ALLOCATE( rmesh(Nr) )
  ALLOCATE( drmesh(Nr) )
  CALL gen_rmesh_exp(r_min, r_max, a, Nr, rmesh)
  CALL gen_drmesh_exp(r_min, r_max, a, Nr, drmesh)

  ! Setup potentials
  ALLOCATE( V(Nr) )
  ALLOCATE( Vmid(Nr-1) )

  Z = 1
  n = 1
  l = 0

  V(:) = -Z/rmesh(:)

  CALL get_midpoints( Nr, rmesh, V, Vmid )

  ALLOCATE( u1(Nr) )
  ALLOCATE( u2(Nr) )

  Emin = -5000.d0
  Emax = 0.d0
  E = -1000.d0

  ctp = Nr ! no perturbation correction

  DO iterShoot = 1,50

    dE = abs(Emax - Emin)
    WRITE(*,*)
    WRITE(*,*) 'iterShoot = ', iterShoot
    WRITE(*,*) 'E = ', E
    WRITE(*,*)
    IF( dE < 1d-10 ) THEN
      WRITE(*,*) 'Converged !!!'
      EXIT
    ENDIF 

    CALL integrate_rschroed_rk4( Nr, l, Z, E, rmesh, V, Vmid, u1, u2, imax )
    !CALL integrate_rproblem_outward( Nr, l, E, rmesh, drmesh, V, &
    !    Z, 0.d0, .false., u1, u2, imax )

    WRITE(*,*) 'imax = ', imax

    !CALL get_min_idx( imax, u1(1:imax), minidx )
    !WRITE(*,*) 'minidx = ', minidx

    !IF( minidx <= 0 ) THEN 
    !  WRITE(1003,*) u1
    !  WRITE(1004,*) u2
    !  WRITE(*,*) "ERROR: The wavefunction doesn't have a peak"
    !  WRITE(*,*) 'u1 and u2 are written anyway'
    !  STOP
    !ENDIF 
  
    ! Trim the wavefunction after the last minimum:
    !u1(minidx:) = 0.d0
    !u2(minidx:) = 0.d0
  
    CALL get_n_nodes( imax, u1(1:imax), nnodes )
    WRITE(*,*) 'DEBUG: nnodes = ', nnodes

    ! use bisection:
    IF( nnodes /= n-l-1 .or. ctp == Nr .or. imax < ctp ) THEN 

      WRITE(*,'(1x,A,2F18.10)') 'DEBUG: before, Emin, Emax = ', Emin, Emax

      CALL is_E_above(n, l, nnodes, isbig)
      !isbig = nnodes > n-l-1
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
    ENDIF 

  ENDDO 


    ! Normalize the wavefunction: (non-relativistic case)
    !CALL integrate_trapz_7(Nr, drmesh, u1**2, S)
    !S = sqrt(abs(S))
    !WRITE(*,*) 'normalization = ', S

    !IF( S > 0.d0 ) THEN 
    !  u1(:) = u1(:) / S
    !  u2(:) = u2(:) / S
    !ELSE 
      ! This would happen if the function is zero, but we already check this
      ! above (converged == 4), so we fail laudly here.
      !  CALL stop_error("solve_radial_eigenproblem: zero function")
    !ENDIF

    !WRITE(1003,*) u1
    !WRITE(1004,*) u2

    !WRITE(*,*) 'Pass here'

  DEALLOCATE( u1 )
  DEALLOCATE( u2 )
  DEALLOCATE( rmesh )
  DEALLOCATE( drmesh )
  DEALLOCATE( V )
  DEALLOCATE( Vmid )

END PROGRAM 

