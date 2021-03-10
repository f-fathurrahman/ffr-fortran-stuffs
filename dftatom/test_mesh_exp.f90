program test_mesh_exp
  use mesh, only: mesh_exp
  implicit none
  real(8) :: r_min, r_max, a
  integer :: N
  real(8), allocatable :: r(:)
  integer :: i
  
  N = 10
  r_min = 0.d0
  r_max = 50.d0
  a = 1d9
  allocate(r(N+1))
  
  r = mesh_exp(r_min, r_max, a, N)

  do i = 1,N+1
    write(*,'(1x,I8,ES18.10)') i, r(i)
  enddo

  deallocate(r)
end