!=========================================================================
subroutine init_gaussian_general(nx,ny,nz,alpha,x0,ga)
  use m_definitions
  use m_tools
  use m_gaussian
  implicit none
  integer,intent(in) :: nx,ny,nz
  real(dp),intent(in) :: alpha,x0(3)
  type(gaussian),intent(out) :: ga
  !=====
  real(dp) :: double_factorial

 ga%nx = nx
 ga%ny = ny
 ga%nz = nz
 ga%am = nx + ny + nz
 ga%amc = orbital_momentum_name(ga%am)

 ga%alpha = alpha

 ga%common_norm_factor = ( 2.0_dp / pi )**0.75_dp &
                 * 2.0_dp**ga%am * ga%alpha**( 0.25_dp * ( 2.0_dp*ga%am + 3.0_dp ) ) 
                 
 ga%norm_factor = ga%common_norm_factor &
                   / SQRT( REAL( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) , dp ) )

 ga%x0(:) = x0(:)

end subroutine init_gaussian_general