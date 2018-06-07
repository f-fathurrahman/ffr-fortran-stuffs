!=========================================================================
subroutine numerical_nucleus(ga,gb)
  use m_definitions
  use m_tools,only: coeffs_gausslegint
  use m_gaussian
 implicit none
 type(gaussian),intent(in) :: ga,gb
!=====
 integer,parameter  :: nx=200
 integer,parameter  :: n1=86
 real(dp),parameter :: rmax=5.
 real(dp)           :: dx,rtmp,x(3)
 integer            :: ix,iy,iz
 integer            :: iangular,na
 real(dp)           :: wu(nx),u(nx),weight
 real(dp)           :: x1(n1),y1(n1),z1(n1),w1(n1)
!=====
  real(dp) :: eval_gaussian

 dx = rmax/REAL(nx,dp)

 write(stdout,*) 'hydrogen atom in zero'

 !
 ! spherical integration
 ! radial part with Gauss-Legendre
 call coeffs_gausslegint(-1.0_dp,1.0_dp,u,wu,nx)
 !
 ! Transformation from [-1;1] to [0;+\infty[
 ! taken from M. Krack JCP 1998
 wu(:) = wu(:) * ( 1.0_dp / log(2.0_dp) / ( 1.0_dp - u(:) ) )
 u(:)  = log( 2.0_dp / (1.0_dp - u(:) ) ) / log(2.0_dp)
! call ld0038(x1,y1,z1,w1,na)
 call ld0086(x1,y1,z1,w1,na)

 rtmp=0.0_dp
 do ix=1,nx
   do iangular=1,n1
     x(1) = u(ix) * x1(iangular)
     x(2) = u(ix) * y1(iangular)
     x(3) = u(ix) * z1(iangular)
     weight = wu(ix) * w1(iangular) * u(ix)**2 * 4.0_dp * pi

     if( SUM(x(:)**2) < 1.d-8 ) cycle ! skip the singularity

     rtmp = rtmp + eval_gaussian(ga,x) * eval_gaussian(gb,x) * weight  / SQRT(SUM(x(:)**2)) * (-1.0_dp)

   enddo
 enddo


 write(stdout,*) 'check V_ab',rtmp

end subroutine numerical_nucleus