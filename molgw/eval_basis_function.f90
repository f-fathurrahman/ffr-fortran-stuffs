!=========================================================================
function eval_basis_function(bf,x)
  use m_definitions, only : dp
  use m_basis_set
  use m_gaussian
  implicit none
  type(basis_function), intent(in) :: bf
  real(dp), intent(in) :: x(3)
  !
  real(dp) :: eval_basis_function
  !=====
  integer                         :: ig
  !=====

  eval_basis_function=0.0_dp
  do ig=1,bf%ngaussian
    eval_basis_function = eval_basis_function + eval_gaussian(bf%g(ig),x) * bf%coeff(ig)
  enddo

end function eval_basis_function
