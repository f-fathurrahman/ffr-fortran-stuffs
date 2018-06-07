!=========================================================================
subroutine evaluate_wfn_r(nspin,nstate,basis,c_matrix,istate1,istate2,ispin,rr,wfn_i)
 use m_definitions
 use m_basis_set
 implicit none
 integer,intent(in)         :: nspin
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nstate
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 integer,intent(in)         :: istate1,istate2,ispin
 real(dp),intent(in)        :: rr(3)
 real(dp),intent(out)       :: wfn_i(istate1:istate2)
!=====
 real(dp)                   :: basis_function_r(basis%nbf)
!=====

 ! First precalculate all the needed basis function evaluations at point rr
 call calculate_basis_functions_r(basis,rr,basis_function_r)

 ! Then rotate
 wfn_i(istate1:istate2) = MATMUL( basis_function_r(:) , c_matrix(:,istate1:istate2,ispin) )


end subroutine evaluate_wfn_r