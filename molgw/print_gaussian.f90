!=========================================================================
subroutine print_gaussian(ga)
 implicit none
 type(gaussian),intent(in) :: ga
!=====

 write(stdout,*)
 write(stdout,*) '*** output information of a cartesian gaussian function ***'
 write(stdout,'(a30,5x,a1)')      'orbital momentum',ga%amc
 write(stdout,'(a30,2x,3(1x,i3))') 'momentum composition',ga%nx,ga%ny,ga%nz
 write(stdout,'(a30,(1x,f12.6))')  'alpha coeff',ga%alpha
 write(stdout,'(a30,3(1x,f12.6))') 'center of the function',ga%x0(:)
 write(stdout,*)

end subroutine print_gaussian