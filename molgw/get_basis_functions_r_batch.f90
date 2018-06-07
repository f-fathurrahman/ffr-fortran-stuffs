!=========================================================================
subroutine get_basis_functions_r_batch(basis,igrid,nr,basis_function_r)
  use m_definitions
  use m_basis_set
  use m_dft_grid
 implicit none
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid,nr
 real(dp),intent(out)       :: basis_function_r(basis%nbf,nr)
!=====
!=====

 ! Check if the batch had been fully precalculated
 ! else calculate it now.
 if( igrid <= ngrid_stored .AND. igrid+nr-1 <= ngrid_stored ) then
   basis_function_r(:,:) = bfr(:,igrid:igrid+nr-1)
 else
   call calculate_basis_functions_r_batch(basis,nr,rr_grid(:,igrid:igrid+nr-1),basis_function_r)
 endif

end subroutine get_basis_functions_r_batch