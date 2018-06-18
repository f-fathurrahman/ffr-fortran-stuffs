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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the density given the Kohn-Sham orbitals.
subroutine build_rho(w, rho)
  use states
  use mesh

  implicit none
  real(8), intent(in)  :: w(N, N, N_wf)
  real(8), intent(out) :: rho(N, N)

  integer :: ix, iy, i

  rho = 0.0_8
  do i = 1, N_occ
     do ix = 1, N
        do iy = 1, N
           rho(ix, iy) = rho(ix, iy) + 2.0_8*w(ix, iy, i)**2
        enddo
     enddo
  enddo

end subroutine build_rho

subroutine zbuild_rho(w, rho)
  use states
  use mesh

  implicit none
  complex(8), intent(in) :: w(N, N, N_wf)
  real(8), intent(out) :: rho(N, N)

  integer :: ix, iy, i

  rho = 0.0_8
  do i = 1, N_occ
     do ix = 1, N
        do iy = 1, N
           rho(ix, iy) = rho(ix, iy) + 2.0_8*conjg(w(ix, iy, i))*w(ix, iy, i)
        enddo
     enddo
  enddo

end subroutine zbuild_rho
