!=========================================================================
subroutine setup_exchange_longrange(print_matrix_,nbf,p_matrix,exchange_ij,eexchange)
  use m_definitions
 use m_eri
  use m_timing
  use m_inputparam, only : nspin, spin_fact
 implicit none
 logical,intent(in)   :: print_matrix_
 integer,intent(in)   :: nbf
 real(dp),intent(in)  :: p_matrix(nbf,nbf,nspin)
 real(dp),intent(out) :: exchange_ij(nbf,nbf,nspin)
 real(dp),intent(out) :: eexchange
!=====
 integer              :: ibf,jbf,kbf,lbf,ispin
!=====
  logical :: negligible_basispair
  real(dp) :: eri_lr

 write(stdout,*) 'Calculate Long-Range Exchange term'
 call start_clock(timing_exchange)


 exchange_ij(:,:,:)=0.0_dp

 do ispin=1,nspin
   do jbf=1,nbf
     do ibf=1,nbf
       do lbf=1,nbf
         do kbf=1,nbf
           !
           ! symmetry (ik|lj) = (ki|lj) has been used to loop in the fast order
           exchange_ij(ibf,jbf,ispin) = exchange_ij(ibf,jbf,ispin) &
                      - eri_lr(kbf,ibf,lbf,jbf) * p_matrix(kbf,lbf,ispin) / spin_fact 
         enddo
       enddo
     enddo
   enddo
 enddo


 eexchange = 0.5_dp*SUM(exchange_ij(:,:,:)*p_matrix(:,:,:))

 call stop_clock(timing_exchange)

end subroutine setup_exchange_longrange
