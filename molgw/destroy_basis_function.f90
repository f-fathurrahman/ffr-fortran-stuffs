!=========================================================================
subroutine destroy_basis_function(bf)
  use m_basis_set
 implicit none
 type(basis_function),intent(inout) :: bf
!=====
 
 deallocate(bf%g,bf%coeff)

end subroutine destroy_basis_function