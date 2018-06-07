!=========================================================================
subroutine destroy_basis_set(basis)
  use m_basis_set
 implicit none

 type(basis_set),intent(inout) :: basis
!=====
 integer :: ibf,ishell
!=====

! do ibf=1,basis%nbf_cart
!   call destroy_basis_function(basis%bfc(ibf))
! enddo
! do ibf=1,basis%nbf
!   call destroy_basis_function(basis%bff(ibf))
! enddo
 deallocate(basis%bfc)
 deallocate(basis%bff)
 do ishell=1,basis%nshell
   if(ALLOCATED(basis%shell(ishell)%alpha)) deallocate( basis%shell(ishell)%alpha )
   if(ALLOCATED(basis%shell(ishell)%coeff)) deallocate( basis%shell(ishell)%coeff )
 enddo
 deallocate(basis%shell)

end subroutine destroy_basis_set
