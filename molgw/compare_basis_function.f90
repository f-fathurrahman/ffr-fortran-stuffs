!=========================================================================
function compare_basis_function(bf1,bf2) result(same_basis_function)
  use m_basis_set
  implicit none
  logical :: same_basis_function
  type(basis_function),intent(in) :: bf1,bf2
  !=====
  integer :: ig
  !=====
   logical :: compare_gaussian

  same_basis_function = .TRUE.

  ! DO NOT compare the following commented fields. Not really necessary...
  ! bf1%shell_index
  ! bf1%amc
  if( bf1%am            /= bf2%am                        ) same_basis_function = .FALSE.
  if( bf1%nx            /= bf2%nx                        ) same_basis_function = .FALSE.
  if( bf1%ny            /= bf2%ny                        ) same_basis_function = .FALSE.
  if( bf1%nz            /= bf2%nz                        ) same_basis_function = .FALSE.
  if( bf1%iatom         /= bf2%iatom                     ) same_basis_function = .FALSE.
  if( ANY(ABS(bf1%x0(:) - bf2%x0(:)) > 1.0e-5_dp )       ) same_basis_function = .FALSE.
  if( bf1%ngaussian     /= bf2%ngaussian                 ) same_basis_function = .FALSE.

  ! If the basis functions already differs, then exit immediately
  if( .NOT. same_basis_function ) return

  do ig=1,bf1%ngaussian
    same_basis_function = same_basis_function .AND. compare_gaussian(bf1%g(ig),bf2%g(ig))
  enddo
  if( ANY(ABS(bf1%coeff(:) - bf2%coeff(:)) > 1.0e-5_dp ) ) same_basis_function = .FALSE.
 
end function compare_basis_function