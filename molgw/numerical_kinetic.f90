!=========================================================================
subroutine numerical_kinetic(ga,gb)
 implicit none
 type(gaussian),intent(in) :: ga,gb
!=====
 integer,parameter  :: nx=100
 real(dp),parameter :: rmax=10.0_dp
 real(dp),parameter :: dh=1.d-4
 real(dp)           :: dx,rtmp,x(3),dhx(3),dhy(3),dhz(3)
 integer            :: ix,iy,iz
!=====

 dx = rmax/REAL(nx,dp)
 dhx(:) = 0.0_dp
 dhx(1) = dh
 dhy(:) = 0.0_dp
 dhy(2) = dh
 dhz(:) = 0.0_dp
 dhz(3) = dh

 rtmp=0.0_dp
 do ix=1,nx
   x(1) = ( REAL(ix,dp)/REAL(nx,dp) - 0.5 ) * rmax
   do iy=1,nx
     x(2) = ( REAL(iy,dp)/REAL(nx,dp) - 0.5 ) * rmax
     do iz=1,nx
       x(3) = ( REAL(iz,dp)/REAL(nx,dp) - 0.5 ) * rmax
  
       rtmp = rtmp - 0.5 * eval_gaussian(ga,x) * dx**3 &
                    * ( eval_gaussian(gb,x+dhx) + eval_gaussian(gb,x-dhx) &
                       +eval_gaussian(gb,x+dhy) + eval_gaussian(gb,x-dhy) &
                       +eval_gaussian(gb,x+dhz) + eval_gaussian(gb,x-dhz) &
                       - 6.0 * eval_gaussian(gb,x) ) / dh**2 
  
     enddo
   enddo
 enddo

 write(stdout,*) 'check K_ab',rtmp

end subroutine numerical_kinetic