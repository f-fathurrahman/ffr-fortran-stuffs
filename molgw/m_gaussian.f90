!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform operations on the single Cartesian Gaussian functions
!
!=========================================================================
module m_gaussian
 use m_definitions
 use m_mpi
 use m_gos
 use m_cart_to_pure

 ! type containing all the information for one unnormalized cartesian gaussian
 ! x**nx * y**ny * z**nz * exp( - alpha * ( x**2 + y**2 + z**2 ) )
 type gaussian
   integer          :: am
   character(len=1) :: amc
   integer          :: nx,ny,nz
   real(dp)         :: alpha
   real(dp)         :: x0(3)                ! center of the gaussian
   real(dp)         :: norm_factor          ! normalization factor for the gaussian squared
   real(dp)         :: common_norm_factor   ! normalization factor for the gaussian squared
                                            ! without the nx,ny,nz dependence
 end type gaussian

 interface
   subroutine boys_function_c(fnt,n,t) BIND(C,NAME='boys_function_c')
     import :: C_INT,C_DOUBLE
     integer(C_INT),value       :: n
     real(C_DOUBLE),value       :: t
     real(C_DOUBLE),intent(out) :: fnt(0:n)
   end subroutine boys_function_c
 end interface

end module m_gaussian

