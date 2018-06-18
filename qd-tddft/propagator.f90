!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA

subroutine propagator(v, dt)
  use states
  use mesh
  use expo
  implicit none

  real(8), intent(in) :: v(N, N)
  real(8), intent(in) :: dt

  integer :: j
  real(8) :: ex, ec
  real(8), allocatable :: vext(:, :), vx(:, :), vc(:, :), rho(:, :), vh(:, :), vp(:, :)
  complex(8), allocatable :: zwfsp(:, :, :)


  allocate(zwfsp(N, N, N_wf), vp(N, N), vext(N, N), vc(N, N), vh(N, N), vx(N, N), rho(N, N))

  ! This is the first, auxiliary propagation, to get v(t+dt)
  zwfsp = zwfs
  do j = 1, N_occ
     call exponential(v, zwfsp(:, :, j), dt)
  enddo
  call external_pot(vext)
  call zbuild_rho(zwfsp, rho)
  call interaction_pot(rho, vh, vx, vc, ex, ec)
  vp = vext + vh + vx + vc

  do j = 1, N_occ
     call exponential(v, zwfs(:, :, j), dt/2.0_8)
     call exponential(vp, zwfs(:, :, j), dt/2.0_8)
  enddo

  deallocate(zwfsp, vp, vext, vc, rho)
end subroutine propagator
