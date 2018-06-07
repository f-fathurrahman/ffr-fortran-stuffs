!=========================================================================
subroutine kinetic_recurrence(ga,gb,k_ab)
  use m_definitions
  use m_gaussian
  implicit none
  type(gaussian),intent(in)     :: ga,gb
  real(dp),intent(out)          :: k_ab
!=====
 real(dp)             :: zeta_ab,ksi_ab,ab2,fact
 real(dp)             :: p(3),ap(3),bp(3)
 real(dp)             :: s_tmp_x(0:ga%nx,0:gb%nx)
 real(dp)             :: s_tmp_y(0:ga%ny,0:gb%ny)
 real(dp)             :: s_tmp_z(0:ga%nz,0:gb%nz)
 real(dp)             :: k_tmp_x(0:ga%nx,0:gb%nx)
 real(dp)             :: k_tmp_y(0:ga%ny,0:gb%ny)
 real(dp)             :: k_tmp_z(0:ga%nz,0:gb%nz)
 integer              :: ix,iy,iz
 integer              :: ixa,iya,iza
 integer              :: ixb,iyb,izb
 integer              :: ixap,iyap,izap
 integer              :: ixbp,iybp,izbp
!=====
 real(dp)             :: s_ab
!=====

 ! Follow the notation of Obara and Saika, JCP  87 3963 (1986)
 zeta_ab = ga%alpha + gb%alpha
 ksi_ab   = ga%alpha * gb%alpha / zeta_ab
 ab2    = SUM( ( ga%x0(:)-gb%x0(:) )**2 )
 p(:)    = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / zeta_ab
 ap(:) =  p(:) - ga%x0(:)
 bp(:) =  p(:) - gb%x0(:)
 fact = 0.5_dp / zeta_ab


 !
 ! direction X
 !
 s_tmp_x(0,0) = (pi/zeta_ab)**1.5_dp * EXP( - ksi_ab * ab2 )
 k_tmp_x(0,0) = ksi_ab * ( 3.0_dp - 2.0_dp * ksi_ab * ab2 ) * s_tmp_x(0,0)

 do ix=1,ga%nx+gb%nx

   do ixa=0,MIN(ga%nx,ix)
     ixb=ix-ixa
     if(ixb>gb%nx) cycle

     if(ixa>0) then
       ixap=ixa-1

       s_tmp_x(ixap+1,ixb) = ap(1) * s_tmp_x(ixap,ixb)
       if(ixap>0)  s_tmp_x(ixap+1,ixb) = s_tmp_x(ixap+1,ixb) + fact * ixap * s_tmp_x(ixap-1,ixb)
       if(ixb>0)   s_tmp_x(ixap+1,ixb) = s_tmp_x(ixap+1,ixb) + fact * ixb  * s_tmp_x(ixap,ixb-1)

       k_tmp_x(ixap+1,ixb) = ap(1) * k_tmp_x(ixap,ixb) + 2.0_dp * ksi_ab * s_tmp_x(ixap+1,ixb)
       if(ixap>0)  k_tmp_x(ixap+1,ixb) = k_tmp_x(ixap+1,ixb) + fact * ixap * k_tmp_x(ixap-1,ixb) &
                    -  ksi_ab / ga%alpha * ixap * s_tmp_x(ixap-1,ixb)
       if(ixb>0)   k_tmp_x(ixap+1,ixb) = k_tmp_x(ixap+1,ixb) + fact * ixb  * k_tmp_x(ixap,ixb-1)

     else
       ixbp=ixb-1

       s_tmp_x(ixa,ixbp+1) = bp(1) * s_tmp_x(ixa,ixbp)
       if(ixbp>0) s_tmp_x(ixa,ixbp+1) = s_tmp_x(ixa,ixbp+1) + fact * ixbp * s_tmp_x(ixa,ixbp-1)
       if(ixa>0)  s_tmp_x(ixa,ixbp+1) = s_tmp_x(ixa,ixbp+1) + fact * ixa  * s_tmp_x(ixa-1,ixbp)

       k_tmp_x(ixa,ixbp+1) = bp(1) * k_tmp_x(ixa,ixbp) +  2.0_dp * ksi_ab * s_tmp_x(ixa,ixbp+1)
       if(ixbp>0) k_tmp_x(ixa,ixbp+1) = k_tmp_x(ixa,ixbp+1) + fact * ixbp * k_tmp_x(ixa,ixbp-1) &
                    -  ksi_ab / gb%alpha * ixbp * s_tmp_x(ixa,ixbp-1)
       if(ixa>0)  k_tmp_x(ixa,ixbp+1) = k_tmp_x(ixa,ixbp+1) + fact * ixa  * k_tmp_x(ixa-1,ixbp)

     endif

   enddo

 enddo

 !
 ! direction Y
 !
 s_tmp_y(0,0) = s_tmp_x(ga%nx,gb%nx)
 k_tmp_y(0,0) = k_tmp_x(ga%nx,gb%nx)

 do iy=1,ga%ny+gb%ny

   do iya=0,MIN(ga%ny,iy)
     iyb=iy-iya
     if(iyb>gb%ny) cycle

     if(iya>0) then
       iyap=iya-1

       s_tmp_y(iyap+1,iyb) = ap(2) * s_tmp_y(iyap,iyb)
       if(iyap>0)  s_tmp_y(iyap+1,iyb) = s_tmp_y(iyap+1,iyb) + fact * iyap * s_tmp_y(iyap-1,iyb)
       if(iyb>0)   s_tmp_y(iyap+1,iyb) = s_tmp_y(iyap+1,iyb) + fact * iyb  * s_tmp_y(iyap,iyb-1)

       k_tmp_y(iyap+1,iyb) = ap(2) * k_tmp_y(iyap,iyb) + 2.0_dp * ksi_ab * s_tmp_y(iyap+1,iyb)
       if(iyap>0)  k_tmp_y(iyap+1,iyb) = k_tmp_y(iyap+1,iyb) + fact * iyap * k_tmp_y(iyap-1,iyb) &
                    -  ksi_ab / ga%alpha * iyap * s_tmp_y(iyap-1,iyb)
       if(iyb>0)   k_tmp_y(iyap+1,iyb) = k_tmp_y(iyap+1,iyb) + fact * iyb  * k_tmp_y(iyap,iyb-1)

     else
       iybp=iyb-1

       s_tmp_y(iya,iybp+1) = bp(2) * s_tmp_y(iya,iybp)
       if(iybp>0) s_tmp_y(iya,iybp+1) = s_tmp_y(iya,iybp+1) + fact * iybp * s_tmp_y(iya,iybp-1)
       if(iya>0)  s_tmp_y(iya,iybp+1) = s_tmp_y(iya,iybp+1) + fact * iya  * s_tmp_y(iya-1,iybp)

       k_tmp_y(iya,iybp+1) = bp(2) * k_tmp_y(iya,iybp) +  2.0_dp * ksi_ab * s_tmp_y(iya,iybp+1)
       if(iybp>0) k_tmp_y(iya,iybp+1) = k_tmp_y(iya,iybp+1) + fact * iybp * k_tmp_y(iya,iybp-1) &
                    -  ksi_ab / gb%alpha * iybp * s_tmp_y(iya,iybp-1)
       if(iya>0)  k_tmp_y(iya,iybp+1) = k_tmp_y(iya,iybp+1) + fact * iya  * k_tmp_y(iya-1,iybp)

     endif

   enddo

 enddo

 !
 ! direction Z
 !
 s_tmp_z(0,0) = s_tmp_y(ga%ny,gb%ny)
 k_tmp_z(0,0) = k_tmp_y(ga%ny,gb%ny)

 do iz=1,ga%nz+gb%nz

   do iza=0,MIN(ga%nz,iz)
     izb=iz-iza
     if(izb>gb%nz) cycle

     if(iza>0) then
       izap=iza-1

       s_tmp_z(izap+1,izb) = ap(3) * s_tmp_z(izap,izb)
       if(izap>0)  s_tmp_z(izap+1,izb) = s_tmp_z(izap+1,izb) + fact * izap * s_tmp_z(izap-1,izb)
       if(izb>0)   s_tmp_z(izap+1,izb) = s_tmp_z(izap+1,izb) + fact * izb  * s_tmp_z(izap,izb-1)

       k_tmp_z(izap+1,izb) = ap(3) * k_tmp_z(izap,izb) + 2.0_dp * ksi_ab * s_tmp_z(izap+1,izb)
       if(izap>0)  k_tmp_z(izap+1,izb) = k_tmp_z(izap+1,izb) + fact * izap * k_tmp_z(izap-1,izb) &
                    -  ksi_ab / ga%alpha * izap * s_tmp_z(izap-1,izb)
       if(izb>0)   k_tmp_z(izap+1,izb) = k_tmp_z(izap+1,izb) + fact * izb  * k_tmp_z(izap,izb-1)

     else
       izbp=izb-1

       s_tmp_z(iza,izbp+1) = bp(3) * s_tmp_z(iza,izbp)
       if(izbp>0) s_tmp_z(iza,izbp+1) = s_tmp_z(iza,izbp+1) + fact * izbp * s_tmp_z(iza,izbp-1)
       if(iza>0)  s_tmp_z(iza,izbp+1) = s_tmp_z(iza,izbp+1) + fact * iza  * s_tmp_z(iza-1,izbp)

       k_tmp_z(iza,izbp+1) = bp(3) * k_tmp_z(iza,izbp) +  2.0_dp * ksi_ab * s_tmp_z(iza,izbp+1)
       if(izbp>0) k_tmp_z(iza,izbp+1) = k_tmp_z(iza,izbp+1) + fact * izbp * k_tmp_z(iza,izbp-1) &
                    -  ksi_ab / gb%alpha * izbp * s_tmp_z(iza,izbp-1)
       if(iza>0)  k_tmp_z(iza,izbp+1) = k_tmp_z(iza,izbp+1) + fact * iza  * k_tmp_z(iza-1,izbp)

     endif

   enddo

 enddo


 !
 ! overlap is a by-product
 !
 s_ab = s_tmp_z(ga%nz,gb%nz) * ga%norm_factor * gb%norm_factor 

 !
 ! final result
 !
 k_ab = k_tmp_z(ga%nz,gb%nz) * ga%norm_factor * gb%norm_factor 


end subroutine kinetic_recurrence