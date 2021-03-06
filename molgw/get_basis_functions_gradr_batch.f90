!=========================================================================
subroutine get_basis_functions_gradr_batch(basis,igrid,nr,basis_function_gradr)
  use m_definitions
  use m_basis_set
  use m_dft_grid
  implicit none
  type(basis_set),intent(in) :: basis
  integer,intent(in)         :: igrid,nr
  real(dp),intent(out)       :: basis_function_gradr(basis%nbf,nr,3)

 ! Check if the batch had been fully precalculated
 ! else calculate it now.
 if( igrid <= ngrid_stored .AND. igrid+nr-1 <= ngrid_stored ) then
   basis_function_gradr(:,:,:) = bfgr(:,igrid:igrid+nr-1,:)
 else
   call calculate_basis_functions_gradr_batch(basis,nr,rr_grid(:,igrid:igrid+nr-1),basis_function_gradr)
 endif

end subroutine get_basis_functions_gradr_batch