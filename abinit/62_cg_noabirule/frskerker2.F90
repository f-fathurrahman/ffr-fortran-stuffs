!{\src2tex{textfont=tt}}
!!****m* ABINIT/frskerker2
!! NAME
!! frskerker2
!!
!! FUNCTION
!! provide the ability to compute the
!! penalty function and its first derivative associated
!! with some residuals and a real space dielectric function
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! NOTES 
!! this is neither a function nor a subroutine. This is a module
!! It is made of two functions and one init subroutine
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module frskerker2 

  use defs_basis
  use defs_abitypes
  use m_profiling_abi
  use interfaces_32_util        ! THIS IS MANDATORY TO CALL dotproduct
  use interfaces_53_spacepar
  use interfaces_56_recipspace  ! THIS IS MANDATORY TO CALL LAPLACIAN

  implicit none

  !! common variables copied from input
  integer,save,private                  :: nfft,nspden,ngfft(18)
  real(dp),save,allocatable,private    :: deltaW(:,:),mat(:,:),rdielng(:)
  real(dp),save,private                :: gprimd(3,3)
  type(dataset_type),pointer,save,private  :: dtset_ptr
  type(MPI_type),pointer,save,private  :: mpi_enreg_ptr
  !! common variables computed
  logical,save,private :: ok=.false.
! *************************************************************************

contains
!!***

!!****f* frskerker2/frskerker2__init
!! NAME
!! frskerker2__init
!!
!! FUNCTION
!! initialisation subroutine
!! Copy every variables required for the energy calculation
!! Allocate the required memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      prcrskerker2
!!
!! CHILDREN
!!      laplacian
!!
!! SOURCE
  subroutine frskerker2__init(dtset_in,mpi_enreg_in,nfft_in,ngfft_in,nspden_in,rdielng_in,deltaW_in,gprimd_in,mat_in )


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'frskerker2__init'
!End of the abilint section

    implicit none
!Arguments ------------------------------------
    type(dataset_type),target,intent(in) :: dtset_in
    integer,intent(in)  :: nfft_in,ngfft_in(18),nspden_in
    real(dp),intent(in) :: deltaW_in(nfft_in,nspden_in),mat_in(nfft_in,nspden_in)
    real(dp),intent(in) :: rdielng_in(nfft_in)
    real(dp),intent(in)  :: gprimd_in(3,3)
    type(MPI_type),target,intent(in)  :: mpi_enreg_in

! *************************************************************************
! !allocation and data transfer
! !Thought it would have been more logical to use the privates intrinsic of the module as
! !input variables it seems that it is not possible...
  if(.not.ok) then
   dtset_ptr => dtset_in
   mpi_enreg_ptr => mpi_enreg_in
   nspden=nspden_in
   ngfft=ngfft_in
   nfft=nfft_in
   ABI_ALLOCATE(deltaW,(size(deltaW_in,1),size(deltaW_in,2)))
   ABI_ALLOCATE(mat,(size(mat_in,1),size(mat_in,2)))
   ABI_ALLOCATE(rdielng,(size(rdielng_in)))
   deltaW=deltaW_in
   rdielng=rdielng_in
   mat=mat_in
   gprimd=gprimd_in
   ok = .true.
  end if
 end subroutine frskerker2__init
!!***

!!****f* frskerker2/frskerker2__end
!! NAME
!! frskerker2__end
!!
!! FUNCTION
!! ending subroutine
!! deallocate memory areas
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      prcrskerker2
!!
!! CHILDREN
!!      laplacian
!!
!! SOURCE
  subroutine frskerker2__end()
    use defs_basis
    use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'frskerker2__end'
!End of the abilint section

    implicit none
!Arguments ------------------------------------

!Local variables-------------------------------

! *************************************************************************
  if(ok) then
!  ! set ok to false which prevent using the pf and dpf
   ok = .false.
!  ! free memory
   ABI_DEALLOCATE(deltaW)
   ABI_DEALLOCATE(mat)
   ABI_DEALLOCATE(rdielng)
  end if
 end subroutine frskerker2__end
!!***

!!****f* frskerker2/frskerker2__newvres2
!! NAME
!! frskerker2__newvres2
!!
!! FUNCTION
!! affectation subroutine
!! do the required renormalisation when providing a new value for
!! the density after application of the gradient
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      laplacian
!!
!! SOURCE
  subroutine frskerker2__newvres2(nv1,nv2,x, grad, vrespc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'frskerker2__newvres2'
!End of the abilint section

    implicit none
!Arguments ------------------------------------

    integer,intent(in)    :: nv1,nv2
    real(dp),intent(in)   :: x
    real(dp),intent(inout)::grad(nv1,nv2)
    real(dp),intent(inout)::vrespc(nv1,nv2)
!Local variables-------------------------------

! *************************************************************************
  grad(:,:)=x*grad(:,:)
  vrespc(:,:)=vrespc(:,:)+grad(:,:)
 end subroutine frskerker2__newvres2
!!***

!!****f* frskerker2/frskerker2__pf
!! NAME
!! frskerker2__pf
!!
!! FUNCTION
!! penalty function associated with the preconditionned residuals    !!
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE
  function frskerker2__pf(nv1,nv2,vrespc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'frskerker2__pf'
 use interfaces_56_recipspace
 use interfaces_62_cg_noabirule
!End of the abilint section

    implicit none
!Arguments ------------------------------------
    integer,intent(in)    :: nv1,nv2
    real(dp),intent(in) ::vrespc(nv1,nv2)
    real(dp)            ::frskerker2__pf
!Local variables-------------------------------

    real(dp)            :: buffer1(nv1,nv2),buffer2(nv1,nv2)
    integer             :: ispden
! *************************************************************************

! Added by YP [HAVE_FORTRAN_INTERFACES]

  if(ok) then
   buffer1=vrespc
   call laplacian(gprimd,mpi_enreg_ptr,nfft,nspden,ngfft,dtset_ptr%paral_kgb,rdfuncr=buffer1,laplacerdfuncr=buffer2)
   do ispden=1,nspden
    buffer2(:,ispden)=(vrespc(:,ispden)-((rdielng(:))**2)*buffer2(:,ispden))  &
&    *half  -  deltaW(:,ispden)
   end do
!  pf_rscgres=dotproduct(vrespc,buffer2)*half-dotproduct(vrespc,deltaW)
   frskerker2__pf=dotproduct(nv1,nv2,vrespc,buffer2)
  else
   frskerker2__pf=zero
  end if
 end function frskerker2__pf
!!***

!!****f* frskerker2/frskerker2__dpf
!! NAME
!! frskerker2__dpf
!!
!! FUNCTION
!! derivative of the penalty function
!! actually not the derivative but something allowing minimization
!! at constant density
!! formula from the work of rackowski,canning and wang
!! H*phi - int(phi**2H d3r)phi
!! that is the simple projection of the change on a direction
!! normal to the density changes
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE
  function frskerker2__dpf(nv1,nv2,vrespc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'frskerker2__dpf'
 use interfaces_56_recipspace
!End of the abilint section

    implicit none
!Arguments ------------------------------------
    integer,intent(in) :: nv1,nv2
    real(dp),intent(in)::vrespc(nv1,nv2)
    real(dp)           :: frskerker2__dpf(nv1,nv2)
!Local variables-------------------------------

    real(dp):: buffer1(nv1,nv2),buffer2(nv1,nv2)
    integer :: ispden

! *************************************************************************

  if(ok) then
   buffer1=vrespc
   call laplacian(gprimd,mpi_enreg_ptr,nfft,nspden,ngfft,dtset_ptr%paral_kgb,rdfuncr=buffer1,laplacerdfuncr=buffer2)
   do ispden=1,nspden
    frskerker2__dpf(:,ispden)= vrespc(:,ispden)-deltaW(:,ispden)-((rdielng(:))**2)*buffer2(:,ispden)
   end do
  else
   frskerker2__dpf = zero
  end if
 end function frskerker2__dpf

end module frskerker2
!!***
