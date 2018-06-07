!=========================================================================
subroutine calculate_basis_functions_r(basis,rr,basis_function_r)
  use m_definitions
  use m_basis_set
  USE m_cart_to_pure
 implicit none

 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: rr(3)
 real(dp),intent(out)       :: basis_function_r(basis%nbf)
!=====
 integer              :: gt
 integer              :: i_cart
 integer              :: ishell,ibf1,ibf2,ibf1_cart
 integer              :: ni_cart,li
 real(dp),allocatable :: basis_function_r_cart(:)
!=====
  integer :: get_gaussian_type_tag
  real(dp) :: eval_basis_function
  integer :: number_basis_function_am

 gt = get_gaussian_type_tag(basis%gaussian_type)

 do ishell=1,basis%nshell
   li        = basis%shell(ishell)%am
   ni_cart   = number_basis_function_am('CART',li)
   ibf1      = basis%shell(ishell)%istart
   ibf1_cart = basis%shell(ishell)%istart_cart
   ibf2      = basis%shell(ishell)%iend


   allocate(basis_function_r_cart(ni_cart))

   do i_cart=1,ni_cart
     basis_function_r_cart(i_cart) = eval_basis_function(basis%bfc(ibf1_cart+i_cart-1),rr)
   enddo
   basis_function_r(ibf1:ibf2) = MATMUL(  basis_function_r_cart(:) , cart_to_pure(li,gt)%matrix(:,:) )
   deallocate(basis_function_r_cart)

 enddo


end subroutine calculate_basis_functions_r
