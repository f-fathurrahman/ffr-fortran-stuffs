!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the matrix to transform CART gaussian to PURE (=spherical) gaussian
!
!=========================================================================
module m_cart_to_pure
  use m_definitions, only : dp

  type transform
    real(dp),allocatable :: matrix(:,:)
  end type

  type(transform), allocatable, protected :: cart_to_pure     (:,:)
  type(transform), allocatable, protected :: cart_to_pure_norm(:,:)

  integer,parameter         :: CARTG=1
  integer,parameter         :: PUREG=2

!=========================================================================
end module m_cart_to_pure
