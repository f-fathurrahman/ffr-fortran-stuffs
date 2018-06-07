!=========================================================================
subroutine write_basis_set(unitfile,basis)
  use m_basis_set
 implicit none

 integer,intent(in)         :: unitfile
 type(basis_set),intent(in) :: basis
!=====
 integer :: ibf
!=====

 write(unitfile)  basis%ammax         
 write(unitfile)  basis%nbf           
 write(unitfile)  basis%nbf_cart      
 write(unitfile)  basis%nshell        
 write(unitfile)  basis%gaussian_type
 do ibf=1,basis%nbf_cart
   call write_basis_function(unitfile,basis%bfc(ibf))
 enddo
 
 end subroutine write_basis_set