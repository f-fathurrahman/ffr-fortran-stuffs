!=================================================================
subroutine calculate_eri_3center_eigen(nbf,nstate,c_matrix,mstate_min,mstate_max,nstate_min,nstate_max)
  USE m_definitions
  USE m_inputparam,only: nspin
  USE m_eri
  USE m_timing
  USE m_mpi
  USE m_eri_ao_mo
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
 integer,intent(in)   :: mstate_min,mstate_max,nstate_min,nstate_max
!=====
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp_l(:,:)
 integer              :: ipair
!=====

 call start_clock(timing_eri_3center_eigen)

 write(stdout,'(/,a)') ' Calculate 3-center integrals on eigenstates'


 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center MO integrals',eri_3center_eigen,1,nauxil_3center,mstate_min,mstate_max,nstate_min,nstate_max,1,nspin)
 eri_3center_eigen(:,:,:,:) = 0.0_dp

 call clean_allocate('TMP 3-center ints',eri_3center_tmp_l,nauxil_3center,nbf)

 do klspin=1,nspin

   do lstate=nstate_min,nstate_max
     if( MODULO( lstate - 1 , nproc_ortho ) /= rank_ortho ) cycle

     eri_3center_tmp_l(:,:) = 0.0_dp

     ! Transformation of the first index
     do ipair=1,npair
       kbf = index_basis(1,ipair)
       lbf = index_basis(2,ipair)
       eri_3center_tmp_l(:,kbf) = eri_3center_tmp_l(:,kbf) &
                                       + c_matrix(lbf,lstate,klspin) * eri_3center(:,ipair)
       if( kbf /= lbf ) &
       eri_3center_tmp_l(:,lbf) = eri_3center_tmp_l(:,lbf) &
                                       + c_matrix(kbf,lstate,klspin) * eri_3center(:,ipair)
     enddo


     ! Transformation of the second index
     eri_3center_eigen(:,mstate_min:mstate_max,lstate,klspin) = MATMUL( eri_3center_tmp_l(:,:) , c_matrix(:,mstate_min:mstate_max,klspin) )

   enddo

 enddo ! klspin

 call clean_deallocate('TMP 3-center ints',eri_3center_tmp_l)

 call xsum_ortho(eri_3center_eigen)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen
