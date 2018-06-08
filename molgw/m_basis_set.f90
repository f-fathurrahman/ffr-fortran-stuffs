!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to deal with the basis set and basis functions
! It works for both wavefunctions basis sets and auxiliary basis sets
!
!=========================================================================
MODULE m_basis_set
  USE m_definitions
  USE m_gaussian

  TYPE basis_function
   ! This basis function belongs to a shell of basis functions
   ! with the same exponents and angular momentum
   INTEGER :: shell_index
   INTEGER :: index_in_shell !
   INTEGER :: am ! Angular momentum number: l=0, 1, 2, 3 ...
   CHARACTER(len=1) :: amc ! Angular momentum letter: s, p, d, f ...
   INTEGER :: nx,ny,nz  ! Angular momentum for cartesian gaussians
   INTEGER :: mm    ! Angular momentum for pure gaussians
   INTEGER :: iatom ! Centered on which atom
   REAL(dp) :: x0(3) ! Coordinates of the gaussian center
   INTEGER :: ngaussian ! Number of primitive gausssians
   TYPE(gaussian), ALLOCATABLE :: g(:) ! The primitive gaussian functions
   REAL(dp), ALLOCATABLE :: coeff(:) ! Their mixing coefficients
 END TYPE

 TYPE shell_type
   INTEGER :: ishell
   INTEGER :: am
   INTEGER :: ng
   REAL(dp), ALLOCATABLE :: alpha(:)
   REAL(dp), ALLOCATABLE :: coeff(:)
   REAL(dp)             :: x0(3)
   INTEGER :: iatom
   ! index of the shell's basis functions in the final basis set
   INTEGER :: istart,iend
   ! index of the shell's basis functions in the cartesian basis set
   INTEGER :: istart_cart,iend_cart
 END TYPE  shell_type



  ! A shell is a group of basis functions sharing: 
  ! the same center, 
  ! the same exponents, 
  ! the same mixing coefficients 
  ! and the same angular momentum

  !
  ! A basis set is a list of basis functions
  TYPE basis_set
    !
    ! The list
    INTEGER :: ammax ! Maximum angular momentum contained in the basis set
    INTEGER :: nbf ! Number of basis functions in the basis set
    INTEGER :: nbf_cart  ! Number of underlying Cartesian functions in the basis set
    INTEGER :: nshell ! Number of shells in the basis sets
    CHARACTER(len=4) :: gaussian_type   ! CART or PURE
    TYPE(basis_function), ALLOCATABLE :: bfc(:) ! Cartesian basis function
    TYPE(basis_function), ALLOCATABLE :: bff(:) ! Final basis function (can be Cartesian or Pure)
    TYPE(shell_type), ALLOCATABLE     :: shell(:)
  END TYPE basis_set

!=========================================================================
END MODULE m_basis_set



