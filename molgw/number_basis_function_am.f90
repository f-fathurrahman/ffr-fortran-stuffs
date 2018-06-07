!=========================================================================
pure function number_basis_function_am(gaussian_type,am)
 implicit none
 character(len=4),intent(in) :: gaussian_type
 integer,intent(in)          :: am
 integer                     :: number_basis_function_am
!=====

 select case(gaussian_type)
 case('CART')
   number_basis_function_am = ( ( am + 1 ) * ( am + 2 ) ) / 2
 case('PURE')
   number_basis_function_am = 2 * am + 1
 end select

end function number_basis_function_am
