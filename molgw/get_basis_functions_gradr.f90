!=========================================================================
subroutine get_basis_functions_gradr(basis,igrid,basis_function_gradr)
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid
 real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
!=====
 real(dp)                   :: rr(3)
!=====

 if( igrid <= ngrid_stored ) then
   basis_function_gradr(:,:) = TRANSPOSE(bfgr(:,igrid,:))
 else
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_gradr(basis,rr,basis_function_gradr)
 endif

end subroutine get_basis_functions_gradr