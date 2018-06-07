!=========================================================================
function eval_gaussian(ga,x)
 implicit none
 type(gaussian),intent(in) :: ga
 real(dp),intent(in) :: x(3)
 real(dp) :: eval_gaussian
!=====
 real(dp) :: dx(3),dx2
!=====

 dx(:) = x(:) - ga%x0(:)
 dx2 = SUM( dx(:)**2 )

 eval_gaussian =  dx(1)**ga%nx &
                * dx(2)**ga%ny &
                * dx(3)**ga%nz &
                * EXP( - ga%alpha * dx2 )

 !
 ! normalize the cartesian gaussian
 eval_gaussian = eval_gaussian * ga%norm_factor

end function eval_gaussian