!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the basic operations related to the quadrature needed to evaluate the
! exchange-correlation integrals
!
!=========================================================================
module m_dft_grid
 use m_definitions
 
 !
 ! Grid definition
 integer :: ngrid
 integer :: nradial
 integer :: nangular_fine
 integer :: nangular_coarse

 real(dp),allocatable :: rr_grid(:,:)
 real(dp),allocatable :: w_grid(:)
 real(dp),allocatable :: bf_rad2(:)

 real(dp),parameter :: pruning_limit = 0.75_dp    ! in terms of covalent radius

 real(dp),parameter :: aa = 0.64_dp ! Scuseria value

 real(dp),parameter :: TOL_WEIGHT = 1.0e-14_dp
 real(dp),parameter :: TOL_BF     = 1.0e-08_dp

 !
 ! Function evaluation storage
 integer :: batch_size_
 integer :: ngrid_stored
 real(dp),allocatable :: bfr(:,:)
 real(dp),allocatable :: bfgr(:,:,:)

!=========================================================================
end module m_dft_grid
!=========================================================================
