!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the method to prepare and store the 2-, 3-, and 4-center Coulomb integrals
!
!=========================================================================
module m_eri
 use m_definitions

 REAL(dp), PARAMETER :: TOO_LOW_EIGENVAL=1.0e-6_dp

 !
 ! max length of a record in the ERI file
 INTEGER, PARAMETER :: line_length=1000

 REAL(dp) :: TOL_INT


 REAL(dp), ALLOCATABLE :: eri_4center(:)
 REAL(dp), ALLOCATABLE :: eri_4center_lr(:)
 REAL(dp), ALLOCATABLE :: eri_3center(:,:)
 REAL(dp), ALLOCATABLE :: eri_3center_lr(:,:)

 LOGICAL, ALLOCATABLE :: negligible_shellpair(:,:)
 INTEGER, ALLOCATABLE :: index_pair_1d(:)
 INTEGER, ALLOCATABLE :: index_basis(:,:)
 INTEGER, ALLOCATABLE :: index_shellpair(:,:)
 INTEGER :: nshellpair

 INTEGER, ALLOCATABLE :: shell_bf(:)

 INTEGER :: nbf_eri         ! local copy of nbf
 INTEGER :: nsize           ! size of the eri_4center array
 INTEGER :: npair         ! number of independent pairs (i,j) with i<=j 

 integer :: nauxil_3center     ! size of the 3-center matrix
                               ! may differ from the total number of 3-center integrals due to
                               ! data distribution
 integer :: nauxil_3center_lr  ! size of the 3-center matrix
                               ! may differ from the total number of 3-center integrals due to
                               ! data distribution

 real(dp), ALLOCATABLE  :: eri_3center_sca(:,:)

! Parallelization information for the auxiliary basis
 integer,ALLOCATABLE :: iproc_ibf_auxil(:)
 integer,ALLOCATABLE :: ibf_auxil_g(:)       ! auxil bf index from local to global
 integer,ALLOCATABLE :: ibf_auxil_l(:)       ! auxil bf index from global to local

! Parallelization information for the auxiliary basis (LR part)
 integer,ALLOCATABLE :: iproc_ibf_auxil_lr(:)
 integer,ALLOCATABLE :: ibf_auxil_g_lr(:)
 integer,ALLOCATABLE :: ibf_auxil_l_lr(:)

!=========================================================================
end module m_eri

