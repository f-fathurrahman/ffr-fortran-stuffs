!=========================================================================
SUBROUTINE calculate_basis_functions_laplr(basis,rr,basis_function_gradr,basis_function_laplr)
  USE m_definitions
  USE m_basis_set
  USE m_cart_to_pure
  implicit none
  type(basis_set),intent(in) :: basis
  real(dp),intent(in)        :: rr(3)
  real(dp),intent(out)       :: basis_function_gradr(3,basis%nbf)
  real(dp),intent(out)       :: basis_function_laplr(3,basis%nbf)
!=====
 integer              :: gt
 integer              :: ishell,ibf1,ibf2,ibf1_cart
 integer              :: i_cart
 integer              :: ni_cart,li
 real(dp),allocatable :: basis_function_gradr_cart(:,:)
 real(dp),allocatable :: basis_function_laplr_cart(:,:)
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
    function eval_basis_function_lapl(bf,x)
      use m_definitions, only : dp
      use m_basis_set
      implicit none
      type(basis_function),intent(in) :: bf
      real(dp),intent(in)             :: x(3)
      real(dp)                        :: eval_basis_function_lapl(3)
    end function
  end interface

 gt = get_gaussian_type_tag(basis%gaussian_type)

 do ishell=1,basis%nshell
   li      = basis%shell(ishell)%am
   ni_cart = number_basis_function_am('CART',li)
   ibf1      = basis%shell(ishell)%istart
   ibf1_cart = basis%shell(ishell)%istart_cart
   ibf2      = basis%shell(ishell)%iend


   allocate(basis_function_gradr_cart(3,ni_cart))
   allocate(basis_function_laplr_cart(3,ni_cart))

   do i_cart=1,ni_cart
     basis_function_gradr_cart(:,i_cart)        = eval_basis_function_grad(basis%bfc(ibf1_cart+i_cart-1),rr)
     basis_function_laplr_cart(:,i_cart)        = eval_basis_function_lapl(basis%bfc(ibf1_cart+i_cart-1),rr)
   enddo

   basis_function_gradr(:,ibf1:ibf2) = MATMUL(  basis_function_gradr_cart(:,:) , cart_to_pure(li,gt)%matrix(:,:) )
   basis_function_laplr(:,ibf1:ibf2) = MATMUL(  basis_function_laplr_cart(:,:) , cart_to_pure(li,gt)%matrix(:,:) )
   deallocate(basis_function_gradr_cart,basis_function_laplr_cart)

 enddo


end subroutine calculate_basis_functions_laplr
