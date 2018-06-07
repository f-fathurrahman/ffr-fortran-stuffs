!=========================================================================
subroutine calculate_basis_functions_r_batch(basis,nr,rr,basis_function_r)
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: nr
 real(dp),intent(in)        :: rr(3,nr)
 real(dp),intent(out)       :: basis_function_r(basis%nbf,nr)
!=====
 integer              :: gt
 integer              :: ir
 integer              :: ishell,ibf1,ibf2,ibf1_cart
 integer              :: i_cart
 integer              :: ni_cart,li
 real(dp),allocatable :: basis_function_r_cart(:,:)
!=====


 gt = get_gaussian_type_tag(basis%gaussian_type)

 do ishell=1,basis%nshell
   li      = basis%shell(ishell)%am
   ni_cart = number_basis_function_am('CART',li)
   ibf1      = basis%shell(ishell)%istart
   ibf1_cart = basis%shell(ishell)%istart_cart
   ibf2      = basis%shell(ishell)%iend

   allocate(basis_function_r_cart(ni_cart,nr))

   do ir=1,nr
     do i_cart=1,ni_cart
       basis_function_r_cart(i_cart,ir) = eval_basis_function(basis%bfc(ibf1_cart+i_cart-1),rr(:,ir))
     enddo
   enddo

   basis_function_r(ibf1:ibf2,:) = MATMUL(  TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)) , basis_function_r_cart(:,:) )
   deallocate(basis_function_r_cart)

 enddo


end subroutine calculate_basis_functions_r_batch