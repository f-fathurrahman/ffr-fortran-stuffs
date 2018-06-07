!=========================================================================
subroutine calculate_basis_functions_gradr_batch(basis,nr,rr,basis_function_gradr)
  use m_definitions
  use m_basis_set
  USE m_cart_to_pure
  implicit none

  type(basis_set),intent(in) :: basis
  integer,intent(in)         :: nr
  real(dp),intent(in)        :: rr(3,nr)
  real(dp),intent(out)       :: basis_function_gradr(basis%nbf,nr,3)
!=====
 integer              :: gt
 integer              :: ir
 integer              :: ishell,ibf1,ibf2,ibf1_cart
 integer              :: i_cart
 integer              :: ni_cart,li
 real(dp),allocatable :: basis_function_gradr_cart(:,:,:)
!=====
  integer :: get_gaussian_type_tag
  integer :: number_basis_function_am
  interface 
    function eval_basis_function_grad(bf,x)
      use m_definitions
      use m_basis_set
      type(basis_function), intent(in) :: bf
      real(dp), intent(in) :: x(3)
      real(dp) :: eval_basis_function_grad(3)
    end function
  end interface

 gt = get_gaussian_type_tag(basis%gaussian_type)

 do ishell=1,basis%nshell
   li      = basis%shell(ishell)%am
   ni_cart = number_basis_function_am('CART',li)
   ibf1      = basis%shell(ishell)%istart
   ibf1_cart = basis%shell(ishell)%istart_cart
   ibf2      = basis%shell(ishell)%iend

   allocate(basis_function_gradr_cart(ni_cart,nr,3))

   do ir=1,nr
     do i_cart=1,ni_cart
       basis_function_gradr_cart(i_cart,ir,:) = eval_basis_function_grad(basis%bfc(ibf1_cart+i_cart-1),rr(:,ir))
     enddo
   enddo

   basis_function_gradr(ibf1:ibf2,:,1) = MATMUL(  TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)) , basis_function_gradr_cart(:,:,1) )
   basis_function_gradr(ibf1:ibf2,:,2) = MATMUL(  TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)) , basis_function_gradr_cart(:,:,2) )
   basis_function_gradr(ibf1:ibf2,:,3) = MATMUL(  TRANSPOSE(cart_to_pure(li,gt)%matrix(:,:)) , basis_function_gradr_cart(:,:,3) )

   deallocate(basis_function_gradr_cart)

 enddo


end subroutine calculate_basis_functions_gradr_batch
