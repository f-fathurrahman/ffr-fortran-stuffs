!=========================================================================
subroutine calculate_eri_2center_scalapack(auxil_basis,rcut)
  use m_definitions
  use m_basis_set
  use m_eri
  use m_timing
  use m_inputparam, only : scalapack_block_min
  use m_scalapack
  use m_eri_calculate
  use m_cart_to_pure
 implicit none
 type(basis_set),intent(in)   :: auxil_basis
 real(dp),intent(in)          :: rcut
!=====
 logical                      :: is_longrange
 integer                      :: ishell,kshell
 integer                      :: n1c,n3c
 integer                      :: ni,nk
 integer                      :: ami,amk
 integer                      :: ibf,kbf
 integer                      :: agt
 integer                      :: info
 integer                      :: ibf_auxil,jbf_auxil
 integer                      :: nauxil_neglect,nauxil_kept
 real(dp)                     :: eigval(auxil_basis%nbf)
 real(dp),allocatable         :: integrals(:,:)
 real(dp)                     :: symmetrization_factor
 real(dp),allocatable         :: eri_2center_sqrt(:,:)
 real(dp),allocatable         :: eri_2center_tmp(:,:)
 integer                      :: mlocal,nlocal
 integer                      :: iglobal,jglobal,ilocal,jlocal
 integer                      :: kglobal,klocal
 integer                      :: desc2center(NDEL)
 logical                      :: skip_shell
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: am1,am3
 integer(C_INT)               :: ng1,ng3
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha3(:)
 real(C_DOUBLE)               :: x01(3),x03(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff3(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====
  integer :: number_basis_function_am
  integer :: get_gaussian_type_tag

 call start_clock(timing_eri_2center)

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut
 agt = get_gaussian_type_tag(auxil_basis%gaussian_type)

 if( .NOT. is_longrange ) then
#ifdef HAVE_SCALAPACK
   write(stdout,'(a,i4,a,i4)') ' 2-center integrals distributed using a SCALAPACK grid (LIBINT): ',nprow_3center,' x ',npcol_3center
#else
   write(stdout,'(a)') ' 2-center integrals (LIBINT)'
#endif
 else
#ifdef HAVE_SCALAPACK
   write(stdout,'(a,i4,a,i4)') ' 2-center LR integrals distributed using a SCALAPACK grid (LIBINT): ',nprow_3center,' x ',npcol_3center
#else
   write(stdout,'(a)') ' 2-center LR integrals (LIBINT)'
#endif
 endif

 !
 ! Enforce that all proc participate to the distribution
 ! Much easier to code then
 if( cntxt_3center <= 0 ) call die('distribution not permitted')



 ! Set mlocal => auxil_basis%nbf
 ! Set nlocal => auxil_basis%nbf
 mlocal = NUMROC(auxil_basis%nbf,block_row,iprow_3center,first_row,nprow_3center)
 nlocal = NUMROC(auxil_basis%nbf,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc2center,auxil_basis%nbf,auxil_basis%nbf,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)

 call clean_allocate('tmp 2-center integrals',eri_2center_tmp,mlocal,nlocal)


 ! Initialization need since we are going to symmetrize the matrix then
 eri_2center_tmp(:,:) = 0.0_dp


 do kshell=1,auxil_basis%nshell
   amk = auxil_basis%shell(kshell)%am
   nk  = number_basis_function_am( auxil_basis%gaussian_type , amk )

   ! Check if this shell is actually needed for the local matrix
   skip_shell = .TRUE.
   do kbf=1,nk
     kglobal = auxil_basis%shell(kshell)%istart + kbf - 1
     skip_shell = skip_shell .AND. .NOT. ( ipcol_3center == INDXG2P(kglobal,block_col,0,first_col,npcol_3center) )
   enddo

   if( skip_shell ) cycle

 
   do ishell=1,auxil_basis%nshell
     ami = auxil_basis%shell(ishell)%am
     ni = number_basis_function_am( auxil_basis%gaussian_type , ami )

     !
     ! Order the angular momenta so that libint is pleased
     ! 1) am3 >= am1
     if( amk < ami ) cycle
     if( amk == ami ) then
       symmetrization_factor = 0.5_dp
     else
       symmetrization_factor = 1.0_dp
     endif

     ! Check if this shell is actually needed for the local matrix
     skip_shell = .TRUE.
     do ibf=1,ni
       iglobal = auxil_basis%shell(ishell)%istart + ibf - 1
       skip_shell = skip_shell .AND. .NOT. ( iprow_3center == INDXG2P(iglobal,block_row,0,first_row,nprow_3center) )
     enddo

     if( skip_shell ) cycle


     am1 = auxil_basis%shell(ishell)%am
     am3 = auxil_basis%shell(kshell)%am
     n1c = number_basis_function_am( 'CART' , ami )
     n3c = number_basis_function_am( 'CART' , amk )
     allocate( int_shell( n1c*n3c ) )
     ng1 = auxil_basis%shell(ishell)%ng
     ng3 = auxil_basis%shell(kshell)%ng
     allocate(alpha1(ng1),alpha3(ng3))
     alpha1(:) = auxil_basis%shell(ishell)%alpha(:) 
     alpha3(:) = auxil_basis%shell(kshell)%alpha(:)
     x01(:) = auxil_basis%shell(ishell)%x0(:)
     x03(:) = auxil_basis%shell(kshell)%x0(:)
     allocate(coeff1(auxil_basis%shell(ishell)%ng))
     allocate(coeff3(auxil_basis%shell(kshell)%ng))
     coeff1(:)=auxil_basis%shell(ishell)%coeff(:) * cart_to_pure_norm(0,agt)%matrix(1,1)
     coeff3(:)=auxil_basis%shell(kshell)%coeff(:) * cart_to_pure_norm(0,agt)%matrix(1,1)


     call libint_2center(am1,ng1,x01,alpha1,coeff1, &
                         am3,ng3,x03,alpha3,coeff3, &
                         rcut_libint,int_shell)


     deallocate(alpha1,alpha3)
     deallocate(coeff1,coeff3)

     call transform_libint_to_molgw(auxil_basis%gaussian_type,ami,amk,int_shell,integrals)

     deallocate(int_shell)

     
     do kbf=1,nk
       kglobal = auxil_basis%shell(kshell)%istart + kbf - 1

       if( ipcol_3center == INDXG2P(kglobal,block_col,0,first_col,npcol_3center) ) then
         klocal = INDXG2L(kglobal,block_col,0,first_col,npcol_3center)
       else
         cycle
       endif

       do ibf=1,ni
         iglobal = auxil_basis%shell(ishell)%istart + ibf - 1

         if( iprow_3center == INDXG2P(iglobal,block_row,0,first_row,nprow_3center) ) then
           ilocal = INDXG2L(iglobal,block_row,0,first_row,nprow_3center)
         else
           cycle
         endif


         eri_2center_tmp(ilocal,klocal) = integrals(ibf,kbf) * symmetrization_factor

       enddo
     enddo
 
     deallocate(integrals)

   enddo   ! ishell
 enddo   ! kshell

 !
 ! Symmetrize and then diagonalize the 2-center integral matrix
 !
#ifdef HAVE_SCALAPACK

 call clean_allocate('2-center integrals sqrt',eri_2center_sqrt,mlocal,nlocal)

 ! B = A
 call PDLACPY('A',auxil_basis%nbf,auxil_basis%nbf,eri_2center_tmp,1,1,desc2center,eri_2center_sqrt,1,1,desc2center)
 ! A = A + B**T
 call PDGEADD('T',auxil_basis%nbf,auxil_basis%nbf,1.0d0,eri_2center_sqrt,1,1,desc2center,1.0d0,eri_2center_tmp,1,1,desc2center)
 ! Diagonalize
 call diagonalize_sca(auxil_basis%nbf,desc2center,eri_2center_tmp,eigval,desc2center,eri_2center_sqrt)
 call clean_deallocate('tmp 2-center integrals',eri_2center_tmp)

#else

 ! Symmetrize
 eri_2center_tmp(:,:) = eri_2center_tmp(:,:) + TRANSPOSE( eri_2center_tmp(:,:) )
 ! Diagonalize
 call diagonalize_scalapack(scalapack_block_min,auxil_basis%nbf,eri_2center_tmp,eigval)
 call move_alloc(eri_2center_tmp,eri_2center_sqrt)

#endif


 !
 ! Skip the too small eigenvalues
 nauxil_kept = COUNT( eigval(:) > TOO_LOW_EIGENVAL )
 nauxil_neglect = auxil_basis%nbf - nauxil_kept
 if( .NOT. is_longrange ) then
   nauxil_2center = nauxil_kept
   ! Prepare the distribution of the 3-center integrals
   ! nauxil_3center variable is now set up
   call distribute_auxil_basis(nauxil_2center)
 else
   nauxil_2center_lr = nauxil_kept
   ! Prepare the distribution of the 3-center integrals
   ! nauxil_3center_lr variable is now set up
   call distribute_auxil_basis_lr(nauxil_2center_lr)
 endif

 !
 ! Now resize the 2-center matrix accordingly
 ! Set mlocal => nauxil_3center
 ! Set nlocal => nauxil_kept < auxil_basis%nbf
 mlocal = NUMROC(auxil_basis%nbf,block_row,iprow_3center,first_row,nprow_3center)
 nlocal = NUMROC(nauxil_kept    ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc_2center,auxil_basis%nbf,nauxil_kept,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)


 if( .NOT. is_longrange ) then
   call clean_allocate('Distributed 2-center integrals',eri_2center,mlocal,nlocal)
 else
   call clean_allocate('Distributed LR 2-center integrals',eri_2center_lr,mlocal,nlocal)
 endif

#ifdef HAVE_SCALAPACK
 call clean_allocate('tmp 2-center integrals',eri_2center_tmp,mlocal,nlocal)
 !
 ! Create a rectangular matrix with only 1 / SQRT( eigval) on a diagonal
 eri_2center_tmp(:,:) = 0.0_dp
 do jlocal=1,nlocal
   jglobal = INDXL2G(jlocal,block_col,ipcol_3center,first_col,npcol_3center)
   do ilocal=1,mlocal
     iglobal = INDXL2G(ilocal,block_row,iprow_3center,first_row,nprow_3center)

     if( iglobal == jglobal + nauxil_neglect ) eri_2center_tmp(ilocal,jlocal) = 1.0_dp / SQRT( eigval(jglobal+nauxil_neglect) )

   enddo
 enddo


 if( .NOT. is_longrange ) then
   call PDGEMM('N','N',auxil_basis%nbf,nauxil_2center,auxil_basis%nbf, &
               1.0_dp,eri_2center_sqrt ,1,1,desc2center,  &
                      eri_2center_tmp,1,1,desc_2center,   &
               0.0_dp,eri_2center    ,1,1,desc_2center)
 else
   call PDGEMM('N','N',auxil_basis%nbf,nauxil_2center_lr,auxil_basis%nbf, &
               1.0_dp,eri_2center_sqrt ,1,1,desc2center,  &
                      eri_2center_tmp,1,1,desc_2center,   &
               0.0_dp,eri_2center_lr ,1,1,desc_2center)
 endif

 call clean_deallocate('tmp 2-center integrals',eri_2center_tmp)

#else

 ilocal = 0
 do jlocal=1,auxil_basis%nbf
   if( eigval(jlocal) < TOO_LOW_EIGENVAL ) cycle
   ilocal = ilocal + 1
   eri_2center_sqrt(:,ilocal) = eri_2center_sqrt(:,jlocal) / SQRT( eigval(jlocal) )
 enddo

 if( .NOT. is_longrange ) then
   do ibf_auxil=1,nauxil_3center
     jbf_auxil = ibf_auxil_g(ibf_auxil)
     eri_2center(:,ibf_auxil) = eri_2center_sqrt(:,jbf_auxil)
   enddo
 else
   do ibf_auxil=1,nauxil_3center_lr
     jbf_auxil = ibf_auxil_g_lr(ibf_auxil)
     eri_2center_lr(:,ibf_auxil) = eri_2center_sqrt(:,jbf_auxil)
   enddo
 endif

#endif

 call clean_deallocate('2-center integrals sqrt',eri_2center_sqrt)



 write(stdout,'(/,1x,a)')      'All 2-center integrals have been calculated, diagonalized and stored'
 write(stdout,'(1x,a,i6)')     'Some have been eliminated due to too large overlap ',nauxil_neglect
 write(stdout,'(1x,a,es16.6)') 'because their eigenvalue was lower than:',TOO_LOW_EIGENVAL


 call stop_clock(timing_eri_2center)

end subroutine calculate_eri_2center_scalapack
