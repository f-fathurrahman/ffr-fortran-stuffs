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
! SUBROUTINE EXTERNAL_POT
! =======================
!
! OUTPUT:
!   v [real(8), dimension(n, n)] : The variable where the external potential is put.
!
! This subroutine places in the "v" argument the external potential that
! defines the quantum dot.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine external_pot(v)
  use mesh

  implicit none

  real(8), intent(out) :: v(N, N)    ! the external potential

  v = 0.0
!!!!!! MISSING CODE 2

!!!!!! END OF MISSING CODE

end subroutine external_pot
