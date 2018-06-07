!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform the Atomic Orbital to Molecular Orbital transform
!
!=========================================================================
module m_eri_ao_mo
 use m_definitions

 real(dp),allocatable :: eri_3center_eigen(:,:,:,:)
 real(dp),allocatable :: eri_3center_eigen_lr(:,:,:,:)
 real(dp),allocatable :: eri_3center_eigen_mixed(:,:,:,:)

!=========================================================================
end module m_eri_ao_mo
