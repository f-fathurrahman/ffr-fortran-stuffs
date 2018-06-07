!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to deal with the basis set and basis functions
! It works for both wavefunctions basis sets and auxiliary basis sets
!
!=========================================================================
module m_basis_set
  use m_definitions
  use m_gaussian

  type basis_function
   ! This basis function belongs to a shell of basis functions
   ! with the same exponents and angular momentum
   integer :: shell_index
   integer :: index_in_shell !
   integer :: am ! Angular momentum number: l=0, 1, 2, 3 ...
   character(len=1) :: amc ! Angular momentum letter: s, p, d, f ...
   integer :: nx,ny,nz  ! Angular momentum for cartesian gaussians
   integer :: mm    ! Angular momentum for pure gaussians
   integer :: iatom ! Centered on which atom
   real(dp) :: x0(3) ! Coordinates of the gaussian center
   integer :: ngaussian ! Number of primitive gausssians
   type(gaussian), allocatable :: g(:) ! The primitive gaussian functions
   real(dp), allocatable :: coeff(:) ! Their mixing coefficients
 end type

 type shell_type
   integer              :: ishell
   integer              :: am
   integer              :: ng
   real(dp),allocatable :: alpha(:)
   real(dp),allocatable :: coeff(:)
   real(dp)             :: x0(3)
   integer :: iatom
  ! index of the shell's basis functions in the final basis set
   integer :: istart,iend
   ! index of the shell's basis functions in the cartesian basis set
   integer :: istart_cart,iend_cart
 end type shell_type


  !
  ! A basis set is a list of basis functions
  type basis_set
    !
    ! The list
    integer :: ammax ! Maximum angular momentum contained in the basis set
    integer :: nbf ! Number of basis functions in the basis set
    integer :: nbf_cart  ! Number of underlying Cartesian functions in the basis set
    integer :: nshell ! Number of shells in the basis sets
    ! A shell is a group of basis functions sharing: 
    ! the same center, 
    ! the same exponents, 
    ! the same mixing coefficients 
    ! and the same angular momentum
    character(len=4)                 :: gaussian_type   ! CART or PURE
    type(basis_function),allocatable :: bfc(:) ! Cartesian basis function
    type(basis_function),allocatable :: bff(:) ! Final basis function (can be Cartesian or Pure)
    type(shell_type),allocatable     :: shell(:)
  end type basis_set

!=========================================================================
end module m_basis_set



