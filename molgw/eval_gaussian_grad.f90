!=========================================================================
function eval_gaussian_grad(ga,x)
 implicit none
 type(gaussian),intent(in) :: ga
 real(dp),intent(in) :: x(3)
 real(dp) :: eval_gaussian_grad(3)
!=====
 real(dp) :: dx(3),dx2
 real(dp) :: g_x,g_y,g_z,gp_x,gp_y,gp_z
!=====

 dx(:) = x(:) - ga%x0(:)
 dx2 = SUM( dx(:)**2 )

 g_x = dx(1)**ga%nx
 g_y = dx(2)**ga%ny
 g_z = dx(3)**ga%nz

 gp_x = -2.0_dp * ga%alpha * dx(1)**(ga%nx+1)
 if(ga%nx>0) gp_x = gp_x + REAL(ga%nx,dp) * dx(1)**(ga%nx-1)
 gp_y = -2.0_dp * ga%alpha * dx(2)**(ga%ny+1)
 if(ga%ny>0) gp_y = gp_y + REAL(ga%ny,dp) * dx(2)**(ga%ny-1)
 gp_z = -2.0_dp * ga%alpha * dx(3)**(ga%nz+1)
 if(ga%nz>0) gp_z = gp_z + REAL(ga%nz,dp) * dx(3)**(ga%nz-1)

 eval_gaussian_grad(1) = gp_x * g_y  * g_z 
 eval_gaussian_grad(2) = g_x  * gp_y * g_z 
 eval_gaussian_grad(3) = g_x  * g_y  * gp_z

 !
 ! multiply by the common exponential factor
 ! and normalize the cartesian gaussian
 eval_gaussian_grad(:) = eval_gaussian_grad(:) * ga%norm_factor * EXP( - ga%alpha * dx2 )


end function eval_gaussian_grad