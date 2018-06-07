!=========================================================================
subroutine prepare_eri(basis)
  use m_definitions
  use m_basis_set
  use m_eri
  use m_inputparam,only: integral_level
  implicit none
  !===== 
  type(basis_set),intent(in) :: basis
  !===== 
  logical            :: file_exists
  !===== 

  nbf_eri = basis%nbf

  select case(integral_level)
  case(low)       ! accuracy not guaranted, just for quick test runs
    TOL_INT = 1.0e-04_dp
  case(medium)    ! 10 meV accuracy on potentials
    TOL_INT = 1.0e-06_dp
  case(high)      !  1 meV accuracy on potentials
    TOL_INT = 1.0e-08_dp
  case(very_high) ! almost perfect potentials
    TOL_INT = 1.0e-10_dp
  case(insane)    ! No screening of any integral
    TOL_INT = 0.0_dp
  case default
    call die('integration quality not recognized')
  end select
  write(stdout,'(/,a,es9.2)') ' Tolerance on integrals set to ',TOL_INT


 if(.NOT.ALLOCATED(negligible_shellpair)) then
   call setup_shell_index(basis)
   allocate(negligible_shellpair(basis%nshell,basis%nshell))
   call identify_negligible_shellpair(basis)
   call setup_shellpair(basis)
   call setup_basispair()
 endif


 nsize = (npair*(npair+1))/2

end subroutine prepare_eri
