!=========================================================================
subroutine level_shifting_up(nbf,nstate,s_matrix,c_matrix,occupation,level_shifting_energy,hamiltonian)
  use m_definitions
  use m_inputparam, only: nspin
 implicit none
 integer,intent(in)     :: nbf,nstate
 real(dp),intent(in)    :: s_matrix(nbf,nbf)
 real(dp),intent(in)    :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)    :: occupation(nstate,nspin)
 real(dp),intent(in)    :: level_shifting_energy
 real(dp),intent(inout) :: hamiltonian(nbf,nbf,nspin)
!=====
 integer  :: ispin
 integer  :: ibf,istate
 real(dp) :: sqrt_level_shifting(nstate)
 real(dp) :: matrix_tmp(nbf,nbf)
!=====

 write(stdout,'(/,a)')     ' Level shifting switched on'
 write(stdout,'(a,f12.6)') '   energy shift (eV):',level_shifting_energy * Ha_eV

 if( level_shifting_energy < 0.0_dp ) then
   call die('level_shifting_energy has to be positive!')
 endif


 do ispin=1,nspin
   !
   ! Shift up empty states only
   do istate=1,nstate
     if( occupation(istate,ispin) < completely_empty ) then
       sqrt_level_shifting(istate) = SQRT( level_shifting_energy )
     else
       sqrt_level_shifting(istate) = 0.0_dp
     endif
   enddo
   forall(istate=1:nstate)
     matrix_tmp(:,istate) =  c_matrix(:,istate,ispin) * sqrt_level_shifting(istate)
   end forall

   ! 
   ! M = C * E * tC
   matrix_tmp(:,:) = MATMUL( matrix_tmp(:,1:nstate) , TRANSPOSE(matrix_tmp(:,1:nstate)) )
   ! M = S * M * S
   matrix_tmp(:,:) = MATMUL( s_matrix , MATMUL( matrix_tmp , s_matrix ) )

   ! Finally update the total hamiltonian
   hamiltonian(:,:,ispin) = hamiltonian(:,:,ispin) + matrix_tmp(:,:)

 enddo
 

end subroutine level_shifting_up
