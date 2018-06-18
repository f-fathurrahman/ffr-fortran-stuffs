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
! MODULE MESH
! ===========
!
! This simple module holds the information that describe the KS states:
! the number of orbitals, and the variables that hold them, used by most
! of the programs and subroutines.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module states

  integer :: N_occ    ! Number of occupied orbitals (occupation number is alwasy two)
  integer :: N_empty  ! Number of unoccupied orbitals.
  integer :: N_wf     ! Total (it should alwasy be N_occ + N_wf)

  ! The next two variables hold the wavefunctions; wfs is real, and is used 
  ! by the "gs" program, whereas zwfs is complex and used by the "td" program.
  real(8), allocatable    :: wfs(:, :, :)
  complex(8), allocatable :: zwfs(:, :, :)

end module states
