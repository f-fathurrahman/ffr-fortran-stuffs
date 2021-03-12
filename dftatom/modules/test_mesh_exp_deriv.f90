program test_mesh_exp_deriv
  use mesh, only: mesh_exp_deriv
  implicit none
  real(8) :: r_min, r_max, a
  integer :: N
  real(8), allocatable :: dr(:)
  integer :: i
  
  N = 11
  r_min = 0.d0
  r_max = 50.d0
  a = 1d9
  allocate(dr(N+1))
  
  dr = mesh_exp_deriv(r_min, r_max, a, N)

  do i = 1,N+1
    write(*,'(1x,I8,ES18.10)') i, dr(i)
  enddo

  deallocate(dr)
end