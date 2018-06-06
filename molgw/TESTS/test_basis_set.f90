PROGRAM test_basis_set

  USE m_atoms, ONLY : natom_basis
  USE m_basis_set

  IMPLICIT NONE 
  CHARACTER(len=4) :: gaussian_type
  CHARACTER(len=100) :: basis_path
  CHARACTER(len=100), ALLOCATABLE :: basis_name(:)
  CHARACTER(len=100), ALLOCATABLE :: ecp_basis_name(:)
  TYPE(basis_set) :: basis

  CALL header()

  CALL do_init_atoms()

  gaussian_type = 'cart'
  basis_path = '/home/efefer/WORKS/MOLGW/molgw-1.F/basis/'

  WRITE(*,*) 'natom_basis = ', natom_basis

  ALLOCATE( basis_name(natom_basis) )
  ALLOCATE( ecp_basis_name(natom_basis) )

  basis_name(:) = '3-21G'
  ecp_basis_name(:) = ''

  CALL init_basis_set(basis_path,basis_name,ecp_basis_name,gaussian_type,basis)

  CALL destroy_atoms()

  DEALLOCATE( basis_name )
  DEALLOCATE( ecp_basis_name )


CONTAINS 

  SUBROUTINE do_init_atoms()

    USE m_atoms
    INTEGER :: natom_read, nghost_read
    REAL(8), ALLOCATABLE :: zatom_read(:), x_read(:,:)
    LOGICAL :: calculate_forces

    natom_read = 1
    nghost_read = 0
    calculate_forces = .false.
    ALLOCATE( zatom_read(natom_read + nghost_read) )
    ALLOCATE( x_read(3,natom_read + nghost_read) )

    zatom_read(1) = 47.d0
    x_read(:,1) = (/ 0.d0, 0.d0, 0.d0 /)
  
    CALL init_atoms( natom_read, nghost_read, zatom_read, x_read, calculate_forces )

    CALL output_positions()
  END SUBROUTINE 



END PROGRAM 
