!=========================================================================
subroutine overlap_three_basis_function(bf1,bf2,bf3,overlap)
  use m_definitions, only : dp
  use m_basis_set
 implicit none
 type(basis_function),intent(in) :: bf1,bf2,bf3
 real(dp),intent(out)            :: overlap
!=====
 type(basis_function)            :: bf12
 integer                         :: ig,jg
 real(dp)                        :: overlap_one_gaussian
!=====

 if(mod(bf1%nx+bf2%nx+bf3%nx,2)==1) then
   overlap=0.0_dp
   return
 endif
 if(mod(bf1%ny+bf2%ny+bf3%ny,2)==1) then
   overlap=0.0_dp
   return
 endif
 if(mod(bf1%nz+bf2%nz+bf3%nz,2)==1) then
   overlap=0.0_dp
   return
 endif
 !
 ! first multiply the two first basis functions
 call basis_function_prod(bf1,bf2,bf12)

 !
 ! then overlap the product and the third basis function
 call overlap_basis_function(bf12,bf3,overlap)

 !
 ! don't forget to destroy it, else memory is leaking
 call destroy_basis_function(bf12)


end subroutine overlap_three_basis_function