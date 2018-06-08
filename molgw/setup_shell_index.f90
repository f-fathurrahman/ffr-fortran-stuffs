!=========================================================================
subroutine setup_shell_index(basis)
  use m_basis_set
  use m_eri
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer :: ibf
!=====

 allocate(shell_bf(basis%nbf))
 shell_bf(:) = basis%bff(:)%shell_index

end subroutine setup_shell_index
