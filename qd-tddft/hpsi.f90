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
! SUBROUTINE HPSI
! ===============
!
! INPUT:
!   v [real(8), dimension(n, n)] : Kohn-Sham potential.
!   f [real(8), dimension(n, n)] : wave function on which the KS Hamiltonian is applied.
! --------
! OUTPUT:
!   hf [real(8), dimension(n, n)] : the resulting wavefunction: |hf> = H|f>
!
! This subroutine operates the Kohn-Sham Hamiltonian (the kinetic term, plus
! the local potential v that should contain the sum of all the terms -- external,
! Hartree and exchange and correlatioin) on a given wavefunction f, and puts
! the result of hf.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hpsi(v, f, hf)
  use mesh

  implicit none

  real(8), intent(in)  :: v(N, N)
  real(8), intent(in)  :: f(N, N)    ! function whose hamiltonian is going to be calculated
  real(8), intent(out) :: hf(N, N) ! the results

  hf = 0.0
!!!!!! MISSING CODE 5

!!!!!! END OF MISSING CODE

end subroutine hpsi




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE ZHPSI
! ===============
!
! INPUT:
!   v [real(8), dimension(n, n)] : Kohn-Sham potential.
!   f [complex(8), dimension(n, n)] : wave function on which the KS Hamiltonian is applied.
! --------
! OUTPUT:
!   hf [complex(8), dimension(n, n)] : the resulting wavefunction: |hf> = H|f>
!
! This subroutine operates the Kohn-Sham Hamiltonian (the kinetic term, plus
! the local potential v that should contain the sum of all the terms -- external,
! Hartree and exchange and correlatioin) on a given wavefunction f, and puts
! the result of hf.
!
! NOTES:
! ------
! It is a repetition of previous subroutine hpsi, but that operates on
! complex wavefunctions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zhpsi(v, f, hf)
  use mesh

  implicit none

  real(8), intent(in)  :: v(N, N)
  complex(8), intent(in)  :: f(N, N)    ! function whose hamiltonian is going to be calculated
  complex(8), intent(inout) :: hf(N, N) ! the results

  hf = (0.0, 0.0)
!!!!!! MISSING CODE 6

!!!!!! END OF MISSING CODE

end subroutine zhpsi
