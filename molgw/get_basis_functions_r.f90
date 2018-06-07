!=========================================================================
subroutine get_basis_functions_r(basis,igrid,basis_function_r)
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: igrid
 real(dp),intent(out)       :: basis_function_r(basis%nbf)
!=====
 real(dp)                   :: rr(3)
!=====

 if( igrid <= ngrid_stored ) then
   basis_function_r(:) = bfr(:,igrid)
 else
   rr(:) = rr_grid(:,igrid)
   call calculate_basis_functions_r(basis,rr,basis_function_r)
 endif

end subroutine get_basis_functions_r