!=========================================================================
subroutine calculate_quadrupole_basis(basis,quadrupole_basis)
  use m_definitions
  use m_basis_set
  use m_cart_to_pure
 implicit none
 type(basis_set),intent(in)         :: basis
 real(dp),allocatable,intent(out)   :: quadrupole_basis(:,:,:,:)
!=====
 integer              :: gt
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2,ibf1_cart,jbf1_cart
 integer              :: li,lj,ni_cart,nj_cart,i_cart,j_cart
 integer              :: idir,jdir
 real(dp),allocatable :: quadrupole_cart(:,:,:,:)
!=====
 integer :: unitfile,var_i
 integer :: get_gaussian_type_tag
 integer :: number_basis_function_am

 gt = get_gaussian_type_tag(basis%gaussian_type)

 allocate(quadrupole_basis(basis%nbf,basis%nbf,3,3))


 do jshell=1,basis%nshell
   lj        = basis%shell(jshell)%am
   nj_cart   = number_basis_function_am('CART',lj)
   jbf1      = basis%shell(jshell)%istart
   jbf1_cart = basis%shell(jshell)%istart_cart
   jbf2      = basis%shell(jshell)%iend

   do ishell=1,basis%nshell
     li        = basis%shell(ishell)%am
     ni_cart   = number_basis_function_am('CART',li)
     ibf1      = basis%shell(ishell)%istart
     ibf1_cart = basis%shell(ishell)%istart_cart
     ibf2      = basis%shell(ishell)%iend


     allocate(quadrupole_cart(ni_cart,nj_cart,3,3))

     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         call basis_function_quadrupole(basis%bfc(ibf1_cart+i_cart-1),basis%bfc(jbf1_cart+j_cart-1),quadrupole_cart(i_cart,j_cart,:,:))
       enddo
     enddo

     do jdir=1,3
       do idir=1,3
         quadrupole_basis(ibf1:ibf2,jbf1:jbf2,idir,jdir) = MATMUL( TRANSPOSE( cart_to_pure(li,gt)%matrix(:,:) ) , &
               MATMUL(  quadrupole_cart(:,:,idir,jdir) , cart_to_pure(lj,gt)%matrix(:,:) ) )
       enddo
     enddo

     deallocate(quadrupole_cart)

   enddo
 enddo


end subroutine calculate_quadrupole_basis
