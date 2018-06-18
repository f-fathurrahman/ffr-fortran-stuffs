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

real(8) function energy(eigenval, rho, vh, vc, vx, ec, ex) result(e)
  use states
  use mesh

  implicit none
  real(8), intent(in) :: eigenval(N_wf)
  real(8), intent(in) :: ec, ex, rho(N, N), vh(N, N), vc(N, N), vx(N, N)

  e = 2.0*sum(eigenval(1:N_occ)) + ec + ex
  e = e - sum(rho(:, :)*(vx(:, :) + vc(:, :) + vh(:, :)/2.0))*delta**2

end function energy
