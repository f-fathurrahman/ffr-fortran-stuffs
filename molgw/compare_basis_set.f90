!=========================================================================
function compare_basis_set(basis1,basis2) result(same_basis_set)
  use m_basis_set
  implicit none

  logical                    :: same_basis_set
  type(basis_set),intent(in) :: basis1,basis2
  !=====
  integer :: ibf
  !=====
  logical :: compare_basis_function

  same_basis_set = .TRUE.
 
  if( basis1%ammax         /= basis2%ammax         )  same_basis_set = .FALSE.
  if( basis1%nbf           /= basis2%nbf           )  same_basis_set = .FALSE.
  if( basis1%nbf_cart      /= basis2%nbf_cart      )  same_basis_set = .FALSE.
  if( basis1%nshell        /= basis2%nshell        )  same_basis_set = .FALSE.
  if( basis1%gaussian_type /= basis2%gaussian_type )  same_basis_set = .FALSE.

  ! If the basis sets already differs, then exit immediately
  if( .NOT. same_basis_set ) return

  do ibf=1,basis1%nbf
    same_basis_set = same_basis_set .AND. compare_basis_function(basis1%bfc(ibf),basis2%bfc(ibf))
  enddo

end function compare_basis_set
