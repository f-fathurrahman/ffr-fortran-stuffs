!=========================================================================
subroutine product_gaussian(ga,gb,gprod)
  USE m_gaussian
 implicit none
 type(gaussian),intent(in) :: ga,gb
 type(gaussian),intent(out) :: gprod
!=====

 if( ANY( ABS(ga%x0(:) - gb%x0(:)) > 1.d-6 ) ) call die('different positions not implemented for product')

 call init_gaussian_general(ga%nx+gb%nx,ga%ny+gb%ny,ga%nz+gb%nz,ga%alpha+gb%alpha,ga%x0,gprod)

end subroutine product_gaussian
