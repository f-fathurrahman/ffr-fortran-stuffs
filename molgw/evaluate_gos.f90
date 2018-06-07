!=========================================================================
subroutine evaluate_gos(ga,gb,qvec,gos_ab)
  USE m_definitions
  USE m_gaussian
  USE m_gos
 implicit none
!=====
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(in)       :: qvec(3)
 complex(dp),intent(out)   :: gos_ab
!=====
 complex(dp) :: sumx,sumy,sumz
 complex(dp) :: fx,gx
 complex(dp) :: fy,gy
 complex(dp) :: fz,gz
 complex(dp) :: factor
 real(dp)    :: aa,bb,ab
 integer     :: ip
!=====

 aa = ga%alpha
 bb = gb%alpha
 ab = aa * bb / ( aa + bb )
 fx = (2.0_dp * aa * bb * ( gb%x0(1) - ga%x0(1) ) + im * aa * qvec(1) ) / (aa + bb)
 gx = (2.0_dp * aa * bb * ( ga%x0(1) - gb%x0(1) ) + im * bb * qvec(1) ) / (aa + bb)
 fy = (2.0_dp * aa * bb * ( gb%x0(2) - ga%x0(2) ) + im * aa * qvec(2) ) / (aa + bb)
 gy = (2.0_dp * aa * bb * ( ga%x0(2) - gb%x0(2) ) + im * bb * qvec(2) ) / (aa + bb)
 fz = (2.0_dp * aa * bb * ( gb%x0(3) - ga%x0(3) ) + im * aa * qvec(3) ) / (aa + bb)
 gz = (2.0_dp * aa * bb * ( ga%x0(3) - gb%x0(3) ) + im * bb * qvec(3) ) / (aa + bb)

 !
 ! x summation
 sumx = 0.0_dp
 do ip = 1, gos(ga%nx,gb%nx)%np
   sumx = sumx + gos(ga%nx,gb%nx)%mu(ip)                      &
                  * fx**gos(ga%nx,gb%nx)%alpha(ip)            &
                  * gx**gos(ga%nx,gb%nx)%beta(ip)             &
                  / (2.0_dp * aa)**gos(ga%nx,gb%nx)%gamma(ip) &
                  / (2.0_dp * bb)**gos(ga%nx,gb%nx)%delta(ip) &
                  * (2.0_dp * ab)**gos(ga%nx,gb%nx)%epsilon(ip)
 enddo

 !
 ! y summation
 sumy = 0.0_dp
 do ip = 1, gos(ga%ny,gb%ny)%np
   sumy = sumy + gos(ga%ny,gb%ny)%mu(ip)                      &
                  * fy**gos(ga%ny,gb%ny)%alpha(ip)            &
                  * gy**gos(ga%ny,gb%ny)%beta(ip)             &
                  / (2.0_dp * aa)**gos(ga%ny,gb%ny)%gamma(ip) &
                  / (2.0_dp * bb)**gos(ga%ny,gb%ny)%delta(ip) &
                  * (2.0_dp * ab)**gos(ga%ny,gb%ny)%epsilon(ip)
 enddo

 !
 ! z summation
 sumz = 0.0_dp
 do ip = 1, gos(ga%nz,gb%nz)%np
   sumz = sumz + gos(ga%nz,gb%nz)%mu(ip)                      &
                  * fz**gos(ga%nz,gb%nz)%alpha(ip)            &
                  * gz**gos(ga%nz,gb%nz)%beta(ip)             &
                  / (2.0_dp * aa)**gos(ga%nz,gb%nz)%gamma(ip) &
                  / (2.0_dp * bb)**gos(ga%nz,gb%nz)%delta(ip) &
                  * (2.0_dp * ab)**gos(ga%nz,gb%nz)%epsilon(ip)
 enddo

 factor = ( pi / ( aa + bb ) )**1.5_dp * EXP( -ab * SUM( (ga%x0(:) - gb%x0(:))**2 ) )  &
           * EXP( ( im * DOT_PRODUCT( qvec(:) , aa * ga%x0(:) + bb * gb%x0(:) )        &
                     - 0.25_dp * SUM(qvec(:)**2) ) / ( aa + bb ) )

 gos_ab = factor * sumx * sumy * sumz * ga%norm_factor * gb%norm_factor


end subroutine evaluate_gos
