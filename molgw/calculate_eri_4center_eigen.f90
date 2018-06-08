!=================================================================
subroutine calculate_eri_4center_eigen(nbf,nstate,c_matrix,istate,ijspin,eri_eigenstate_i)
  use m_definitions
 use m_inputparam, only: nspin
  USE m_eri_ao_mo
  USE m_timing
 implicit none

 integer,intent(in)     :: nbf,nstate
 integer,intent(in)     :: istate,ijspin
 real(dp),intent(in)    :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(inout) :: eri_eigenstate_i(nstate,nstate,nstate,nspin)
!=====
 integer,save         :: istate_previous=0
 integer,save         :: ijspin_previous=0
 integer              :: klspin
 integer              :: ibf,jbf,kbf,lbf
 integer              :: jstate,kstate,lstate
 real(dp)             :: eri_tmp3(nbf,nbf,nbf),eri_tmp2(nbf,nbf,nbf)
!=====
  REAL(dp) :: eri

 ! Check if the calculation can be skipped
 if( istate_previous == istate .AND. ijspin_previous == ijspin .AND. ANY(ABS(eri_eigenstate_i(:,:,:,:))>1.0e-6_dp) ) then
   return
 else
   istate_previous = istate
   ijspin_previous = ijspin
 endif


 call start_clock(timing_basis_transform)

 eri_eigenstate_i(:,:,:,:)=0.0_dp
 eri_tmp2(:,:,:)=0.0_dp
 eri_tmp3(:,:,:)=0.0_dp

 do lbf=1,nbf
   do kbf=1,nbf
     do jbf=1,nbf

       do ibf=1,nbf
         eri_tmp3(jbf,kbf,lbf) = eri_tmp3(jbf,kbf,lbf) + eri(ibf,jbf,kbf,lbf) * c_matrix(ibf,istate,ijspin) 
       enddo


     enddo
   enddo
 enddo

 do lbf=1,nbf
   do kbf=1,nbf

     do jstate=1,nstate
       eri_tmp2(jstate,kbf,lbf) = DOT_PRODUCT( eri_tmp3(:,kbf,lbf) , c_matrix(:,jstate,ijspin) )
     enddo

   enddo
 enddo


  
 do klspin=1,nspin

   do lbf=1,nbf
     do kstate=1,nstate
       do jstate=1,nstate
         eri_tmp3(jstate,kstate,lbf) = DOT_PRODUCT( eri_tmp2(jstate,:,lbf) , c_matrix(:,kstate,klspin) )
       enddo
     enddo
   enddo

   do lstate=1,nstate
     do kstate=1,nstate
       do jstate=1,nstate

         eri_eigenstate_i(jstate,kstate,lstate,klspin) = DOT_PRODUCT( eri_tmp3(jstate,kstate,:) , c_matrix(:,lstate,klspin) )

       enddo
     enddo
   enddo

 enddo !klspin


 call stop_clock(timing_basis_transform)

end subroutine calculate_eri_4center_eigen
