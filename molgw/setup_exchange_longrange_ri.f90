!=========================================================================
subroutine setup_exchange_longrange_ri(nbf,nstate,occupation,c_matrix,p_matrix,exchange_ij,eexchange)
  use m_definitions
 use m_eri
  use m_timing
  use m_inputparam, only : nspin, spin_fact
  use m_mpi
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin,istate,ibf_auxil
 integer              :: index_ij
 integer              :: nocc
 real(dp),allocatable :: tmp(:,:)
 real(dp)             :: eigval(nbf)
 integer              :: ipair
!=====
  logical :: negligible_basispair
  real(dp) :: eri

 write(stdout,*) 'Calculate LR Exchange term with Resolution-of-Identity'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:) = 0.0_dp

 allocate(tmp(nauxil_3center_lr,nbf))

 do ispin=1,nspin

   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty)  cycle
     if( MODULO( istate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     tmp(:,:) = 0.0_dp
     do ipair=1,npair
       ibf = index_basis(1,ipair)
       jbf = index_basis(2,ipair)
       tmp(:,ibf) = tmp(:,ibf) + c_matrix(jbf,istate,ispin) * eri_3center_lr(:,ipair)
       if( ibf /= jbf ) &
            tmp(:,jbf) = tmp(:,jbf) + c_matrix(ibf,istate,ispin) * eri_3center_lr(:,ipair)
     enddo

     !exchange_ij(:,:,ispin) = exchange_ij(:,:,ispin) &
     !                   - MATMUL( TRANSPOSE(tmp(:,:)) , tmp(:,:) ) / spin_fact
     ! C = A^T * A + C
     call DSYRK('L','T',nbf,nauxil_3center_lr,-occupation(istate,ispin)/spin_fact,tmp,nauxil_3center_lr,1.0_dp,exchange_ij(:,:,ispin),nbf)
   enddo

   !
   ! Need to symmetrize exchange_ij
   do ibf=1,nbf
     do jbf=ibf+1,nbf
       exchange_ij(ibf,jbf,ispin) = exchange_ij(jbf,ibf,ispin)
     enddo
   enddo

 enddo
 deallocate(tmp)

 call xsum_world(exchange_ij)

 eexchange = 0.5_dp * SUM( exchange_ij(:,:,:) * p_matrix(:,:,:) )

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange_ri
