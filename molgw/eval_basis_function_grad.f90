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

  eval_basis_function_grad(:)=0.0_dp
  do ig=1,bf%ngaussian
    eval_basis_function_grad(:) = eval_basis_function_grad(:) + eval_gaussian_grad(bf%g(ig),x) *   bf%coeff(ig)
  enddo

end function eval_basis_function_grad