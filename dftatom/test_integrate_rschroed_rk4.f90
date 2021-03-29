PROGRAM test

  IMPLICIT NONE 
  REAL(8) :: r_min, r_max, a
  INTEGER :: Nr
  REAL(8), ALLOCATABLE :: rmesh(:), drmesh(:)
  INTEGER :: i
  REAL(8), ALLOCATABLE :: V(:), Vmid(:)
  REAL(8), ALLOCATABLE :: u1(:), u2(:)
  INTEGER :: Z, l
  INTEGER :: imax
  REAL(8) :: E
  REAL(8) :: S
  INTEGER :: minidx, nnodes

  ! Setup radial mesh
  Nr = 1000
  r_min = 1.d-9
  r_max = 50.d0
  a = 1d9
  !a = 1.d0
  ALLOCATE( rmesh(Nr) )
  ALLOCATE( drmesh(Nr) )
  CALL gen_rmesh_exp(r_min, r_max, a, Nr, rmesh)
  CALL gen_drmesh_exp(r_min, r_max, a, Nr, drmesh)

  ! Setup potentials
  ALLOCATE( V(Nr) )
  ALLOCATE( Vmid(Nr-1) )

  Z = 1
  l = 0

  V(:) = -Z/rmesh(:)

  CALL get_midpoints( Nr, rmesh, V, Vmid )

  WRITE(*,*) 'V = ', V(1:4)
  WRITE(*,*) 'Vmid = ', Vmid(1:3)

  WRITE(1000,*) rmesh
  WRITE(1001,*) V
  WRITE(1002,*) Vmid

  !
  ALLOCATE( u1(Nr) )
  ALLOCATE( u2(Nr) )

  E = -0.4d0
  CALL integrate_rschroed_rk4( Nr, l, Z, E, rmesh, V, Vmid, u1, u2, imax )
  !CALL integrate_rproblem_outward( Nr, l, E, rmesh, drmesh, V, &
  !    Z, 0.d0, .false., u1, u2, imax )

  WRITE(*,*)
  WRITE(*,*) 'After integrate_rproblem_outward'
  WRITE(*,*) 'imax = ', imax
  WRITE(*,*) 'Some u1 and u2'
  DO i = 1,5
    WRITE(*,'(1x,I8,2ES20.10)') i, u1(i), u2(i)
  ENDDO 
  WRITE(*,*) '....'
  DO i = imax-5,imax
    WRITE(*,'(1x,I8,2ES20.10)') i, u1(i), u2(i)
  ENDDO 


  CALL get_min_idx( imax, u1(1:imax), minidx )
  WRITE(*,*) 'minidx = ', minidx

  IF( minidx <= 0 ) THEN 
    WRITE(1003,*) u1
    WRITE(1004,*) u2
    WRITE(*,*) "ERROR: The wavefunction doesn't have a peak"
    WRITE(*,*) 'u1 and u2 are written anyway'
    STOP
  ENDIF 
  
  ! Trim the wavefunction after the last minimum:
  !u1(minidx:) = 0.d0
  !u2(minidx:) = 0.d0
  
  ! To make sure the zeros from above are not counted as nodes, we
  ! substract 1 from minidx here:
  CALL get_n_nodes( minidx-1, u1(1:minidx-1), nnodes )
  WRITE(*,*) 'DEBUG: nnodes = ', nnodes

  ! Normalize the wavefunction: (non-relativistic case)
  CALL integrate_trapz_7(Nr, drmesh, u1**2, S)
  S = sqrt(abs(S))
  WRITE(*,*) 'normalization = ', S

  !IF( S > 0.d0 ) THEN 
  !  u1(:) = u1(:) / S
  !  u2(:) = u2(:) / S
  !ELSE 
    ! This would happen if the function is zero, but we already check this
    ! above (converged == 4), so we fail laudly here.
  !  CALL stop_error("solve_radial_eigenproblem: zero function")
  !ENDIF

  WRITE(1003,*) u1
  WRITE(1004,*) u2

  WRITE(*,*) 'Pass here'

  DEALLOCATE( u1 )
  DEALLOCATE( u2 )
  DEALLOCATE( rmesh )
  DEALLOCATE( drmesh )
  DEALLOCATE( V )
  DEALLOCATE( Vmid )

END PROGRAM 

