!=========================================================================
function eval_gaussian_lapl(ga,x)
 implicit none
 type(gaussian),intent(in) :: ga
 real(dp),intent(in) :: x(3)
 real(dp) :: eval_gaussian_lapl(3)
!=====
 real(dp) :: dx(3),dx2
 real(dp) :: g_x,g_y,g_z,gpp_x,gpp_y,gpp_z
!=====

 dx(:) = x(:) - ga%x0(:)
 dx2 = SUM( dx(:)**2 )

 g_x = dx(1)**ga%nx
 g_y = dx(2)**ga%ny
 g_z = dx(3)**ga%nz

 gpp_x = 4.0_dp * ga%alpha**2 * dx(1)**(ga%nx+2) - 2.0_dp * ga%alpha * REAL(ga%nx+1,dp) * dx(1)**ga%nx
 if(ga%nx>0) gpp_x = gpp_x - 2.0_dp * ga%alpha * REAL(ga%nx,dp) * dx(1)**(ga%nx)
 if(ga%nx>1) gpp_x = gpp_x + 2.0_dp * ga%alpha * REAL(ga%nx*(ga%nx-1),dp) * dx(1)**(ga%nx-1)

 gpp_y = 4.0_dp * ga%alpha**2 * dx(2)**(ga%ny+2) - 2.0_dp * ga%alpha * REAL(ga%ny+1,dp) * dx(2)**ga%ny
 if(ga%ny>0) gpp_y = gpp_y - 2.0_dp * ga%alpha * REAL(ga%ny,dp) * dx(2)**(ga%ny)
 if(ga%ny>1) gpp_y = gpp_y + 2.0_dp * ga%alpha * REAL(ga%ny*(ga%ny-1),dp) * dx(2)**(ga%ny-1)

 gpp_z = 4.0_dp * ga%alpha**2 * dx(3)**(ga%nz+2) - 2.0_dp * ga%alpha * REAL(ga%nz+1,dp) * dx(3)**ga%nz
 if(ga%nz>0) gpp_z = gpp_z - 2.0_dp * ga%alpha * REAL(ga%nz,dp) * dx(3)**(ga%nz)
 if(ga%nz>1) gpp_z = gpp_z + 2.0_dp * ga%alpha * REAL(ga%nz*(ga%nz-1),dp) * dx(3)**(ga%nz-1)

 eval_gaussian_lapl(1) = gpp_x * g_y   * g_z 
 eval_gaussian_lapl(2) = g_x   * gpp_y * g_z 
 eval_gaussian_lapl(3) = g_x   * g_y   * gpp_z

 !
 ! multiply by the common exponential factor
 ! and normalize the cartesian gaussian
 eval_gaussian_lapl(:) = eval_gaussian_lapl(:) * ga%norm_factor * EXP( - ga%alpha * dx2 )


end function eval_gaussian_lapl