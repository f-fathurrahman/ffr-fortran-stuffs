!=========================================================================
function eval_basis_function_grad(bf,x)
  use m_definitions, only : dp
  use m_basis_set
  implicit none
  type(basis_function), intent(in) :: bf
  real(dp), intent(in) :: x(3)
  real(dp) :: eval_basis_function_grad(3)
  !=====
  integer :: ig
  !=====
  INTERFACE 
    FUNCTION eval_gaussian_grad(ga,x)
      USE m_definitions
      USE m_gaussian
      type(gaussian),intent(in) :: ga
      real(dp),intent(in) :: x(3)
      real(dp) :: eval_gaussian_grad(3)
    END FUNCTION 
  END INTERFACE 

  eval_basis_function_grad(:)=0.0_dp
  do ig=1,bf%ngaussian
    eval_basis_function_grad(:) = eval_basis_function_grad(:) + eval_gaussian_grad(bf%g(ig),x) *   bf%coeff(ig)
  enddo

end function eval_basis_function_grad
