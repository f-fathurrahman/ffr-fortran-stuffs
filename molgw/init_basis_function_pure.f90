!=========================================================================
subroutine init_basis_function_pure(normalized,ng,am,mm,iatom,x0, &
                                    alpha,coeff,shell_index,index_in_shell,bf)
  
  use m_definitions
  use m_basis_set
  use m_tools, only : orbital_momentum_name
  implicit none
  logical,intent(in)               :: normalized
  integer,intent(in)               :: ng,am,mm,shell_index,iatom,index_in_shell
  real(dp),intent(in)              :: x0(3),alpha(ng)
  real(dp),intent(in)              :: coeff(ng)
  type(basis_function),intent(out) :: bf
  !=====
  integer                          :: ig
  real(dp)                         :: overlap
  !=====

  bf%ngaussian = ng
  allocate(bf%g(bf%ngaussian))
  allocate(bf%coeff(bf%ngaussian))
  bf%nx    = -1
  bf%ny    = -1
  bf%nz    = -1
  bf%am    = am
  bf%mm    = mm
  bf%amc   = orbital_momentum_name(bf%am)
  bf%iatom = iatom
  bf%x0(:) = x0(:)
  bf%shell_index = shell_index
  bf%index_in_shell = index_in_shell

! Do not need this

!  ! All the gaussians of the contraction have the same orbital momentum
!  do ig=1,bf%ngaussian
!    call init_gaussian_general(nx,ny,nz,alpha(ig),x0,bf%g(ig))
!    bf%coeff(ig) = coeff(ig)
!  enddo
! 
!  !
!  ! check the normalization if requested
!  if( normalized ) then
!    call overlap_basis_function(bf,bf,overlap)
!    if( ABS(overlap-1.0_dp) > 2.0e-5_dp ) then
! !     write(stdout,*) 'normalization is different from 1.0',overlap
! !     write(stdout,*) bf%nx,bf%ny,bf%nz
! !     write(stdout,*) 'assuming this is a generalized contraction and rescaling coefficients'
!      bf%coeff(:) = coeff(:) / SQRT( overlap )
!    endif
!  endif
 

end subroutine init_basis_function_pure