!=========================================================================
function eval_basis_function_lapl(bf,x)
  use m_definitions, only : dp
  use m_basis_set
  implicit none
  type(basis_function),intent(in) :: bf
  real(dp),intent(in)             :: x(3)
  real(dp)                        :: eval_basis_function_lapl(3)
  !=====
  integer                         :: ig
  !=====
  INTERFACE 
    function eval_gaussian_lapl(ga,x)
      USE m_definitions
      USE m_gaussian
      type(gaussian),intent(in) :: ga
      real(dp),intent(in) :: x(3)
      real(dp) :: eval_gaussian_lapl(3)
    END FUNCTION 
  END INTERFACE 

  eval_basis_function_lapl(:)=0.0_dp
  do ig=1,bf%ngaussian
    eval_basis_function_lapl(:) = eval_basis_function_lapl(:) + eval_gaussian_lapl(bf%g(ig),x) *   bf%coeff(ig)
  enddo

end function eval_basis_function_lapl
