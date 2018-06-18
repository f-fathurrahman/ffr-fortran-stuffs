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




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE PERTURBATION
!
!
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine perturbation(N_wf, wfs)
  use mesh
  implicit none
  integer, intent(in) :: N_wf
  complex(8), intent(inout) :: wfs(N, N, N_wf)

  integer :: j, ix, iy

  do j = 1, N_wf
     do ix = 1, n
        do iy = 1, n
           wfs(ix, iy, j) = exp((0._8,1._8)*0.01_8*x(ix,iy))*wfs(ix, iy, j)
        enddo
     enddo
  enddo

end subroutine perturbation
