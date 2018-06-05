PROGRAM test_atoms
  
  USE m_atoms
  IMPLICIT NONE 
  INTEGER :: natom_read, nghost_read
  REAL(8), ALLOCATABLE :: zatom_read(:), x_read(:,:)
  LOGICAL :: calculate_forces
  INTEGER :: Ncore_states

  natom_read = 1
  nghost_read = 0
  calculate_forces = .false.
  ALLOCATE( zatom_read(natom_read + nghost_read) )
  ALLOCATE( x_read(3,natom_read + nghost_read) )

  zatom_read(1) = 47.d0
  x_read(:,1) = (/ 0.d0, 0.d0, 0.d0 /)
  
  CALL init_atoms( natom_read, nghost_read, zatom_read, x_read, calculate_forces )

  Ncore_states = atoms_core_states()
  WRITE(*,*) 'Ncore_states = ', Ncore_states

  CALL output_positions()

  CALL destroy_atoms()

  WRITE(*,*)
  WRITE(*,*) 'Pass here ...'

END PROGRAM 
