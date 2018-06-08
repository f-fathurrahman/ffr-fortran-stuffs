============================================
subroutine calculate_eri_3center_scalapack(basis,auxil_basis,rcut)
  USE m_definitions
  USE m_basis_set
  use m_eri
  use m_timing
  use m_scalapack
 implicit none
 type(basis_set),intent(in)   :: basis
 type(basis_set),intent(in)   :: auxil_basis
 real(dp),intent(in)          :: rcut
!=====
 logical                      :: is_longrange
 integer                      :: agt
 integer                      :: ishell,kshell,lshell
 integer                      :: klshellpair
 integer                      :: n1c,n3c,n4c
 integer                      :: ni,nk,nl
 integer                      :: ami,amk,aml
 integer                      :: ibf,kbf,lbf
 integer                      :: ibf_auxil,jbf_auxil
 integer                      :: info
 real(dp),allocatable         :: integrals(:,:,:)
 real(dp),allocatable         :: eri_3center_tmp(:,:)
 integer                      :: klpair_global
 integer                      :: ibf_auxil_local,ibf_auxil_global
 integer                      :: mlocal,nlocal
 integer                      :: iglobal,jglobal,ilocal,jlocal
 integer                      :: kglobal,klocal
 integer                      :: desc3center(NDEL)
 integer                      :: desc3tmp(NDEL)
 integer                      :: desc3final(NDEL)
 integer                      :: nauxil_kept
 logical                      :: skip_shell
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: am1,am3,am4
 integer(C_INT)               :: ng1,ng3,ng4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====
  integer :: index_pair
  integer :: number_basis_function_am
  integer :: get_gaussian_type_tag

 call start_clock(timing_eri_3center)

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut
 agt = get_gaussian_type_tag(auxil_basis%gaussian_type)

 if( .NOT. is_longrange ) then
   nauxil_kept = nauxil_2center
 else
   nauxil_kept = nauxil_2center_lr
 endif

 if( .NOT. is_longrange ) then
   write(stdout,'(/,a)')    ' Calculate and store all the 3-center Electron Repulsion Integrals (LIBINT 3center)'
 else
   write(stdout,'(/,a)')    ' Calculate and store all the LR 3-center Electron Repulsion Integrals (LIBINT 3center)'
 endif
#ifdef HAVE_SCALAPACK
 write(stdout,'(a,i4,a,i4)') ' 3-center integrals distributed using a SCALAPACK grid: ',nprow_3center,' x ',npcol_3center
#endif

 !
 ! Enforce that all proc participate to the distribution
 ! Much easier to code then
 if( cntxt_3center <= 0 ) call die('distribution not permitted')


 ! Set mlocal => auxil_basis%nbf
 ! Set nlocal => npair
 mlocal = NUMROC(auxil_basis%nbf,block_row,iprow_3center,first_row,nprow_3center)
 nlocal = NUMROC(npair          ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3center,auxil_basis%nbf,npair,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)

 !  Allocate the 3-center integral array
 !
 ! 3-CENTER INTEGRALS 
 !
 call clean_allocate('TMP 3-center integrals',eri_3center_tmp,mlocal,nlocal)



 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am
   nk = number_basis_function_am( basis%gaussian_type , amk )
   nl = number_basis_function_am( basis%gaussian_type , aml )


   ! Check if this shell is actually needed for the local matrix
   skip_shell = .TRUE.
   do lbf=1,nl
     do kbf=1,nk
       klpair_global = index_pair(basis%shell(kshell)%istart+kbf-1,basis%shell(lshell)%istart+lbf-1)

       skip_shell = skip_shell .AND. .NOT. ( ipcol_3center == INDXG2P(klpair_global,block_col,0,first_col,npcol_3center) )
     enddo
   enddo

   if( skip_shell ) cycle


   do ishell=1,auxil_basis%nshell

     ami = auxil_basis%shell(ishell)%am
     ni = number_basis_function_am( auxil_basis%gaussian_type , ami )

     ! Check if this shell is actually needed for the local matrix
     skip_shell = .TRUE.
     do ibf=1,ni
       iglobal = auxil_basis%shell(ishell)%istart + ibf - 1
       skip_shell = skip_shell .AND. .NOT. ( iprow_3center == INDXG2P(iglobal,block_row,0,first_row,nprow_3center) )
     enddo

     if( skip_shell ) cycle


     am1 = ami
     am3 = amk
     am4 = aml
     n1c = number_basis_function_am( 'CART' , ami )
     n3c = number_basis_function_am( 'CART' , amk )
     n4c = number_basis_function_am( 'CART' , aml )
     ng1 = auxil_basis%shell(ishell)%ng
     ng3 = basis%shell(kshell)%ng
     ng4 = basis%shell(lshell)%ng
     allocate(alpha1(ng1),alpha3(ng3),alpha4(ng4))
     allocate(coeff1(ng1),coeff3(ng3),coeff4(ng4))
     alpha1(:) = auxil_basis%shell(ishell)%alpha(:) 
     alpha3(:) = basis%shell(kshell)%alpha(:)
     alpha4(:) = basis%shell(lshell)%alpha(:)
     coeff1(:) = auxil_basis%shell(ishell)%coeff(:) * cart_to_pure_norm(0,agt)%matrix(1,1)
     coeff3(:) = basis%shell(kshell)%coeff(:)
     coeff4(:) = basis%shell(lshell)%coeff(:)
     x01(:) = auxil_basis%shell(ishell)%x0(:)
     x03(:) = basis%shell(kshell)%x0(:)
     x04(:) = basis%shell(lshell)%x0(:)


     allocate( int_shell(n1c*n3c*n4c) )

     call libint_3center(am1,ng1,x01,alpha1,coeff1, &
                         am3,ng3,x03,alpha3,coeff3, &
                         am4,ng4,x04,alpha4,coeff4, &
                         rcut_libint,int_shell)

     call transform_libint_to_molgw(auxil_basis%gaussian_type,ami,basis%gaussian_type,amk,aml,int_shell,integrals)

     
     do lbf=1,nl
       do kbf=1,nk
         klpair_global = index_pair(basis%shell(kshell)%istart+kbf-1,basis%shell(lshell)%istart+lbf-1)
         if( ipcol_3center /= INDXG2P(klpair_global,block_col,0,first_col,npcol_3center) ) cycle
         jlocal = INDXG2L(klpair_global,block_col,0,first_col,npcol_3center)

         do ibf=1,ni
           ibf_auxil_global = auxil_basis%shell(ishell)%istart+ibf-1
           if( iprow_3center /= INDXG2P(ibf_auxil_global,block_row,0,first_row,nprow_3center) ) cycle
           ilocal = INDXG2L(ibf_auxil_global,block_row,0,first_row,nprow_3center)

           eri_3center_tmp(ilocal,jlocal) = integrals(ibf,kbf,lbf)

         enddo
       enddo
     enddo


     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,alpha3,alpha4)
     deallocate(coeff1,coeff3,coeff4)

   enddo ! ishell

 enddo ! klshellpair


 ! Set mlocal => nauxil_kept = nauxil_2center OR nauxil_2center_lr
 ! Set nlocal => npair
 mlocal = NUMROC(nauxil_kept,block_row,iprow_3center,first_row,nprow_3center)
 nlocal = NUMROC(npair      ,block_col,ipcol_3center,first_col,npcol_3center)
 call DESCINIT(desc3tmp,nauxil_kept,npair,block_row,block_col,first_row,first_col,cntxt_3center,MAX(1,mlocal),info)




#ifdef HAVE_SCALAPACK
 call clean_allocate('3-center integrals SCALAPACK',eri_3center_sca,mlocal,nlocal)

 if( .NOT. is_longrange ) then
   call PDGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center    ,1,1,desc_2center,  &
                      eri_3center_tmp,1,1,desc3center,   &
               0.0_dp,eri_3center_sca,1,1,desc3tmp)
 else
   call PDGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center_lr ,1,1,desc_2center,  &
                      eri_3center_tmp,1,1,desc3center,   &
               0.0_dp,eri_3center_sca,1,1,desc3tmp)
 endif

 write(stdout,'(a,i8,a,i4)') ' Final 3-center integrals distributed using a SCALAPACK grid: ',nprow_auxil,' x ',npcol_auxil

 if( cntxt_auxil > 0 ) then
   mlocal = NUMROC(nauxil_kept,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   nlocal = NUMROC(npair      ,NBLOCK_AUXIL,ipcol_auxil,first_col,npcol_auxil)
 else
   mlocal = -1
   nlocal = -1
 endif
 call xmax_ortho(mlocal)
 call xmax_ortho(nlocal)

 call clean_deallocate('TMP 3-center integrals',eri_3center_tmp)

 !
 ! Final distribution on Full-Range or Long-Range arrays
 !
 call DESCINIT(desc3final,nauxil_kept,npair,MBLOCK_AUXIL,NBLOCK_AUXIL,first_row,first_col,cntxt_auxil,MAX(1,mlocal),info)
 if( .NOT. is_longrange ) then
   call clean_deallocate('Distributed 2-center integrals',eri_2center)
   call clean_allocate('3-center integrals',eri_3center,mlocal,nlocal)
   call PDGEMR2D(nauxil_kept,npair,eri_3center_sca,1,1,desc3tmp,eri_3center,1,1,desc3final,cntxt_3center)
   !
   ! Propagate to the ortho MPI direction
   if( cntxt_auxil <= 0 ) then
     eri_3center(:,:) = 0.0_dp
   endif
   call xsum_ortho(eri_3center)

   ! Development version
#ifndef SCASCA
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
#endif

 else
   call clean_deallocate('Distributed LR 2-center integrals',eri_2center_lr)
   call clean_allocate('LR 3-center integrals',eri_3center_lr,mlocal,nlocal)
   call PDGEMR2D(nauxil_kept,npair,eri_3center_sca,1,1,desc3tmp,eri_3center_lr,1,1,desc3final,cntxt_3center)
   !
   ! Propagate to the ortho MPI direction
   if( cntxt_auxil <= 0 ) then
     eri_3center_lr(:,:) = 0.0_dp
   endif
   call xsum_ortho(eri_3center_lr)

  ! Development version
#ifndef SCASCA
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
#endif
 endif


#else

 if( .NOT. is_longrange ) then

   call clean_allocate('3-center integrals',eri_3center,mlocal,nlocal)

   call DGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center,auxil_basis%nbf,eri_3center_tmp,auxil_basis%nbf,0.0_dp,eri_3center,nauxil_kept)

   write(stdout,'(a)') ' Final 3-center integrals not distributed! since there is no SCALAPACK'

   call clean_deallocate('Distributed 2-center integrals',eri_2center)

 else

   call clean_allocate('LR 3-center integrals',eri_3center_lr,mlocal,nlocal)

   call DGEMM('T','N',nauxil_kept,npair,auxil_basis%nbf, &
               1.0_dp,eri_2center_lr,auxil_basis%nbf,eri_3center_tmp,auxil_basis%nbf,0.0_dp,eri_3center_lr,nauxil_kept)

   write(stdout,'(a)') ' Final LR 3-center integrals not distributed! since there is no SCALAPACK'

   call clean_deallocate('Distributed 2-center integrals',eri_2center_lr)

 endif
 call clean_deallocate('TMP 3-center integrals',eri_3center_tmp)

 
#endif



 write(stdout,'(/,1x,a)') 'All 3-center integrals have been calculated and stored'


 call stop_clock(timing_eri_3center)


end subroutine calculate_eri_3center_scalapack
