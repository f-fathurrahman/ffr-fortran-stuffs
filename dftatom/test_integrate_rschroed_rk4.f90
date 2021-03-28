PROGRAM test

  IMPLICIT NONE 
  REAL(8) :: r_min, r_max, a
  INTEGER :: Nr
  REAL(8), ALLOCATABLE :: rmesh(:)
  INTEGER :: i
  REAL(8), ALLOCATABLE :: V(:), Vmid(:)
  REAL(8), ALLOCATABLE :: u1(:), u2(:)
  INTEGER :: Z, l
  INTEGER :: imax
  REAL(8) :: E

  ! Setup radial mesh
  Nr = 1000
  r_min = 1.d-9
  r_max = 50.d0
  a = 1d9
  !a = 1.d0
  ALLOCATE( rmesh(Nr) )
  CALL gen_rmesh_exp(r_min, r_max, a, Nr, rmesh)

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

  E = -0.1d0
  CALL integrate_rschroed_rk4( Nr, l, Z, E, rmesh, V, Vmid, u1, u2, imax )

  WRITE(1003,*) u1
  WRITE(1004,*) u2

  WRITE(*,*) 'umax = ', imax
  WRITE(*,*) 'Pass here'

  DEALLOCATE( u1 )
  DEALLOCATE( u2 )
  DEALLOCATE( rmesh )
  DEALLOCATE( V )
  DEALLOCATE( Vmid )

END PROGRAM 

