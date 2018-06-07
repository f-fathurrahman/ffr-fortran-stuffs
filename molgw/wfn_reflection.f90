!=========================================================================
function wfn_reflection(nstate,basis,c_matrix,istate,ispin)
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
 integer                    :: wfn_reflection
!=====
 real(dp) :: phi_tmp1(1),phi_tmp2(1)
 real(dp) :: xtmp1(3),xtmp2(3)
 real(dp) :: proj
!=====

 xtmp1(1) = x(1,1) +  2.0_dp
 xtmp1(2) = x(2,1) +  1.0_dp
 xtmp1(3) = x(3,1) +  3.0_dp
 call evaluate_wfn_r(nspin,nstate,basis,c_matrix,istate,istate,ispin,xtmp1,phi_tmp1)

 proj = DOT_PRODUCT( xtmp1 , xnormal )
 xtmp2(:) = xtmp1(:) -  2.0_dp * proj * xnormal(:)
 call evaluate_wfn_r(nspin,nstate,basis,c_matrix,istate,istate,ispin,xtmp2,phi_tmp2)

 if( ABS(phi_tmp1(1) - phi_tmp2(1))/ABS(phi_tmp1(1)) < 1.0e-6_dp ) then
   wfn_reflection = 1
 else if( ABS(phi_tmp1(1) + phi_tmp2(1))/ABS(phi_tmp1(1)) < 1.0e-6_dp ) then
   wfn_reflection = -1
 else 
   wfn_reflection = 0
 endif


end function wfn_reflection
