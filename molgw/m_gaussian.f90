!=========================================================================
! This file is part of MOLGW.
! Author: Fabien Bruneval
!
! This module contains
! the methods to perform operations on the single Cartesian Gaussian functions
!
!=========================================================================
MODULE m_gaussian
 USE m_definitions

 ! type containing all the information for one unnormalized cartesian gaussian
 ! x**nx * y**ny * z**nz * exp( - alpha * ( x**2 + y**2 + z**2 ) )
 TYPE gaussian
   INTEGER :: am
   CHARACTER(len=1) :: amc
   INTEGER  :: nx,ny,nz
   REAL(dp) :: alpha
   REAL(dp) :: x0(3)                ! center of the gaussian
   REAL(dp) :: norm_factor          ! normalization factor for the gaussian squared
   ! normalization factor for the gaussian squared without the nx,ny,nz dependence
   REAL(dp) :: common_norm_factor
 END TYPE gaussian

 INTERFACE 
   SUBROUTINE boys_function_c(fnt,n,t) BIND(C,NAME='boys_function_c')
     import :: C_INT,C_DOUBLE
     INTEGER(C_INT), VALUE       :: n
     REAL(C_DOUBLE), VALUE       :: t
     REAL(C_DOUBLE), INTENT(out) :: fnt(0:n)
   END SUBROUTINE boys_function_c
 END INTERFACE

END MODULE m_gaussian

