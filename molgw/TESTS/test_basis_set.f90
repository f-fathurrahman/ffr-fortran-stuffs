PROGRAM test_basis_set

  USE m_atoms, ONLY : natom_basis
  USE m_atoms, ONLY : destroy_atoms
  USE m_basis_set

  IMPLICIT NONE 
  CHARACTER(len=4) :: gaussian_type
  CHARACTER(len=100) :: basis_path
  CHARACTER(len=100), ALLOCATABLE :: basis_name(:)
  CHARACTER(len=100), ALLOCATABLE :: ecp_basis_name(:)
  TYPE(basis_set) :: basis
  INTEGER :: ibf

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

  CALL setup_cart_to_pure_transforms()

  WRITE(*,*)
  WRITE(*,*) "gaussian_type = ", basis%gaussian_type
  WRITE(*,*) "nbf = ", basis%nbf
  WRITE(*,*) "nbf_cart = ", basis%nbf_cart
  WRITE(*,*)
  !DO ibf = 1,basis%nbf
  !  CALL print_basis_function(basis%bff(ibf))
  !ENDDO 

  !DO ibf = 1,basis%nbf_cart
  !  CALL print_basis_function(basis%bfc(ibf))
  !ENDDO 

  CALL test_integrals()

  CALL destroy_atoms()

  DEALLOCATE( basis_name )
  DEALLOCATE( ecp_basis_name )


CONTAINS 


  SUBROUTINE test_integrals()
    IMPLICIT NONE 
    REAL(8) :: ovl
    CALL print_basis_function( basis%bfc(1) )
    CALL overlap_basis_function( basis%bfc(1), basis%bfc(1), ovl )
    WRITE(*,*) 'ovl = ', ovl
  END SUBROUTINE 

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
