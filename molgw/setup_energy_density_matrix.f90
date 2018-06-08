!=========================================================================
subroutine setup_energy_density_matrix(nbf,nstate,c_matrix,occupation,energy,q_matrix)
  use m_definitions
  use m_timing
  use m_inputparam, only: nspin
 implicit none
 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)  :: occupation(nstate,nspin)
 real(dp),intent(in)  :: energy(nstate,nspin)
 real(dp),intent(out) :: q_matrix(nbf,nbf)
!=====
 integer :: ispin,ibf,jbf
 integer :: istate
!=====

 call start_clock(timing_density_matrix)
 write(stdout,'(1x,a)') 'Build energy-density matrix'

 q_matrix(:,:) = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) cycle
     call DSYR('L',nbf,occupation(istate,ispin)*energy(istate,ispin),c_matrix(:,istate,ispin),1,q_matrix(:,:),nbf)
   enddo
 enddo


 ! Symmetrize
 do jbf=1,nbf
   do ibf=jbf+1,nbf
     q_matrix(jbf,ibf) = q_matrix(ibf,jbf)
   enddo
 enddo
 call stop_clock(timing_density_matrix)


end subroutine setup_energy_density_matrix
