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

subroutine propagate()
  use states
  use mesh
  implicit none

  real(8), external :: energy
  complex(8), allocatable :: hwf(:, :)
  real(8), allocatable :: rho(:, :)
  real(8), allocatable :: vext(:,:), vh(:,:), vx(:,:), vc(:,:), vtot(:, :) ! the potentials
  real(8), allocatable :: eigenval(:)

  integer :: niter, i, j, ix, iy, dipole_unit
  real(8) :: dt, t, ex, ec, etot, dipole(1:2), prop_time

  allocate(rho(N, N), vext(N, N), vh(N, N), vx(N, N), vc(N, N), vtot(N, N), hwf(N, N))
  allocate(eigenval(N_wf))

  prop_time = 1000.0_8
  dt = 1.0_8
  niter = nint(prop_time/dt)

  dipole_unit = 12
  open(unit = dipole_unit, file="dipole", action ="write", status="replace", form ="formatted")
  close(unit = dipole_unit)
  
  write(*,'(/,a,/)') 'Time dependent propagation follows.'
  do i = 0, niter, 1

     call external_pot(vext)
     call zbuild_rho(zwfs, rho)
     call interaction_pot(rho, vh, vx, vc, ex, ec)
     vtot = vext + vh + vx + vc
     do j = 1, N_wf
        call zhpsi(vtot, zwfs(:, :, j), hwf(:, :))
        eigenval(j) = zdotproduct(zwfs(:, :, j), hwf)
     enddo
        
     etot =  energy(eigenval, rho, vh, vc, vx, ec, ex)

     dipole = 0.0_8
     do ix = 1, n
        do iy = 1, n
           dipole(1) = dipole(1) + x(ix, iy)*rho(ix, iy)
           dipole(2) = dipole(2) + y(ix, iy)*rho(ix, iy)
        enddo
     enddo

     t = i*dt

     open(unit = dipole_unit, file="dipole", position="append")
     write(dipole_unit, '(3e18.8)') t, dipole(1), dipole(2)
     close(unit = dipole_unit)

     write(*, *) i, t, etot

     call propagator(vtot, dt)

  enddo

  close(unit=dipole_unit)
  
  deallocate(rho, vext, vh, vx, vc, vtot, hwf, eigenval)

end subroutine propagate
