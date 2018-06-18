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
! SUBROUTINE INTERACTION_POT
! ==========================
!
! INPUT:
!   rho [real(8), dimension(n, n)] : electronic density.
! ---------
! OUTPUT:
!   vx [real(8), dimension(n, n)] : Hartree potential.
!   vx [real(8), dimension(n, n)] : exchange-energy potential.
!   vc [real(8), dimension(n, n)] : correlation-energy potential.
!   ex [real(8)] : exchange energy.
!   ec [real(8)] : correlation energy.
!
! Given an input density, this subroutine must provide the Hartree, exchange
! and correlation potential, as well as the exchange energy and potential.
! For that purpose, it must call "poisson_solve" and "vxc_lda"
!
! It is useful if the code is built in such a way that is easy to choose
! between full use of LDA, use of LDA for exchange only (null correlation),
! use only of the Hartree term (exchange and correlation are null), and
! the independent particle approximation (everything is null).
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interaction_pot(rho, vh, vx, vc, ex, ec)
  use mesh
  use poisson
  implicit none
  real(8), intent(in) :: rho(n, n)
  real(8), intent(out) :: vh(n, n), vx(n, n), vc(n, n), ex, ec


  vh = 0.0; vx = 0.0; vc = 0.0; ex = 0.0; ec = 0.0
!!!!!! MISSING CODE 4

!!!!!! END OF MISSING CODE

end subroutine interaction_pot





