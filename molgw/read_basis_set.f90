!=========================================================================
subroutine read_basis_set(unitfile,basis)
  use m_basis_set
 implicit none

 integer,intent(in)          :: unitfile
 type(basis_set),intent(out) :: basis
!=====
 integer :: ibf
!=====

 read(unitfile)  basis%ammax
 read(unitfile)  basis%nbf
 read(unitfile)  basis%nbf_cart
 read(unitfile)  basis%nshell
 read(unitfile)  basis%gaussian_type
 allocate(basis%bfc(basis%nbf_cart))
 do ibf=1,basis%nbf_cart
   call read_basis_function(unitfile,basis%bfc(ibf))
 enddo


end subroutine read_basis_set