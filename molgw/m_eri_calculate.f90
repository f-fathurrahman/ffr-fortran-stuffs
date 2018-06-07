!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the method to calculate the 2-, 3-, and 4-center Coulomb integrals
!
!=========================================================================
module m_eri_calculate
 use m_definitions
 use m_scalapack, only : NDEL

 integer :: nauxil_2center     ! size of the 2-center matrix
 integer :: nauxil_2center_lr  ! size of the 2-center LR matrix

 real(dp),allocatable :: eri_2center(:,:)
 real(dp),allocatable :: eri_2center_lr(:,:)
 integer :: desc_2center(NDEL)

!=========================================================================
end module m_eri_calculate
