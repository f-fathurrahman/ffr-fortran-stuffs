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
! PROGRAM TEST_LAPLACIAN
! ======================
!
! This program is built with the purpose of testing the implementation of the
! Laplacian subroutines in the mesh module. It builds a Gaussian distribution
! (normalized so that it integrates to one), and calculates its Laplacian.
! Then it compares it with the exact analytical vales, and outputs the error,
! defined to be <(f_approx - f_exact)|(f_aaprox - f_exact)>
!
! You may try different function, or vary the "hardness" of the function by
! changing the alpha parameter. You may also check how the error depends
! on the order of the discretization of the Laplacian. 
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_laplacian
  use mesh
  implicit none

  integer :: ix, iy
  real(8), allocatable :: rho(:, :), exl(:, :), apl(:, :)
  real(8) :: r2, alpha, z, ex, err
  real(8), parameter :: pi = 3.141592653589793_8

  allocate(rho(n, n), exl(n, n), apl(n, n))


  ! Define the problem density, which is just a Gaussian.
  alpha = 3.0d0
  do ix = 1, n
     do iy = 1, n
        r2 = x(ix, iy)**2 + y(ix, iy)**2
        rho(ix, iy) = exp(-r2/alpha**2)
     enddo
  enddo
  rho(:, :) = rho(:, :)/(alpha**2*pi)

  ! The exact value of the Laplacian of the problem density is put in exl variable
  do ix = 1, n
     do iy = 1, n
        r2 = x(ix, iy)**2 + y(ix, iy)**2
        exl(ix, iy) = (4.0d0/alpha**2)*(r2/alpha**2-1.0d0)*rho(ix, iy)
     enddo
  enddo

  ! Calculate the Laplacian of the problem density throught the laplacian
  ! subroutine, and put the result into apl variable.
  call laplacian(rho, apl)

  ! For visualization purposes:
  call output(rho,'rho')
  call output(exl,'exact_laplacian')
  call output(apl,'approximated_laplacian')

  ! Outputs an estimation of the error of the laplacian subroutine.
  write(*, '(a,es20.8)') 'Error: ', dotproduct(apl-exl,apl-exl)

  deallocate(rho, exl, apl)
end subroutine test_laplacian
