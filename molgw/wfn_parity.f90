!=========================================================================
function wfn_parity(nstate,basis,c_matrix,istate,ispin)
 use m_definitions
 use m_mpi
 use m_atoms
 use m_basis_set
 use m_inputparam
 implicit none
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 integer,intent(in)         :: istate,ispin
 integer                    :: wfn_parity
!=====
 real(dp) :: phi_tmp1(1),phi_tmp2(1),xtmp(3)
!=====

 xtmp(1) = xcenter(1) +  2.0_dp
 xtmp(2) = xcenter(2) +  1.0_dp
 xtmp(3) = xcenter(3) +  3.0_dp
 call evaluate_wfn_r(nspin,nstate,basis,c_matrix,istate,istate,ispin,xtmp,phi_tmp1)
 xtmp(1) = xcenter(1) -  2.0_dp
 xtmp(2) = xcenter(2) -  1.0_dp
 xtmp(3) = xcenter(3) -  3.0_dp
 call evaluate_wfn_r(nspin,nstate,basis,c_matrix,istate,istate,ispin,xtmp,phi_tmp2)

 if( ABS(phi_tmp1(1) - phi_tmp2(1))/ABS(phi_tmp1(1)) < 1.0e-6_dp ) then
   wfn_parity = 1
 else
   wfn_parity = -1
 endif
 

end function wfn_parity