!=================================================================
subroutine calculate_eri_3center_eigen_lr(nbf,nstate,c_matrix)
  USE m_definitions
  USE m_inputparam,only: nspin
  USE m_eri
  USE m_timing
  USE m_mpi
  USE m_eri_ao_mo
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
!=====
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp_l(:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate LR 3-center integrals on eigenstates'


 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('LR 3-center MO integrals',eri_3center_eigen_lr,nauxil_3center_lr,nstate,nstate,nspin)
 eri_3center_eigen_lr(:,:,:,:) = 0.0_dp

 allocate(eri_3center_tmp_l(nauxil_3center_lr,nbf))

 do klspin=1,nspin

   do lstate=1,nstate
     if( MODULO( lstate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     eri_3center_tmp_l(:,:) = 0.0_dp

     ! Transformation of the first index
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       eri_3center_tmp_l(:,kbf) = eri_3center_tmp_l(:,kbf) &
                                       + c_matrix(lbf,lstate,klspin) * eri_3center_lr(:,ipair)
       if( kbf /= lbf )  &
         eri_3center_tmp_l(:,lbf) = eri_3center_tmp_l(:,lbf) &
                                         + c_matrix(kbf,lstate,klspin) * eri_3center_lr(:,ipair)

     enddo

   ! Transformation of the second index
     eri_3center_eigen_lr(:,:,lstate,klspin) = MATMUL( eri_3center_tmp_l(:,:) , c_matrix(:,:,klspin) )

   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp_l)

 call xsum_ortho(eri_3center_eigen)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen_lr
