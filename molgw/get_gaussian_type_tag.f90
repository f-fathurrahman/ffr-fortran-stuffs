!=========================================================================
function get_gaussian_type_tag(gaussian_type)
 implicit none
 character(len=4),intent(in) :: gaussian_type
 integer  :: get_gaussian_type_tag
!=====

 select case(gaussian_type)
 case('CART')
   get_gaussian_type_tag = CARTG
 case('PURE')
   get_gaussian_type_tag = PUREG
 end select
 
end function get_gaussian_type_tag