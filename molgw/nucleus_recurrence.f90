!=========================================================================
subroutine nucleus_recurrence(zatom,c,ga,gb,v_ab)
  use m_definitions
  use m_gaussian
 implicit none
 real(dp),intent(in)       :: zatom,c(3)
 type(gaussian),intent(in) :: ga,gb
 real(dp),intent(out)      :: v_ab
!=====
 real(dp)             :: zeta_ab,ksi_ab,ab2,fact
 real(dp)             :: p(3),ap(3),bp(3),cp(3)
 real(dp)             :: v_tmp_x_m(0:ga%nx,0:gb%nx)
 real(dp)             :: v_tmp_y_m(0:ga%ny,0:gb%ny)
 real(dp)             :: v_tmp_z_m(0:ga%nz,0:gb%nz)
 real(dp)             :: v_tmp_x_mp1(0:ga%nx,0:gb%nx)
 real(dp)             :: v_tmp_y_mp1(0:ga%ny,0:gb%ny)
 real(dp)             :: v_tmp_z_mp1(0:ga%nz,0:gb%nz)
 integer              :: ix,iy,iz
 integer              :: ixa,iya,iza
 integer              :: ixb,iyb,izb
 integer              :: ixap,iyap,izap
 integer              :: ixbp,iybp,izbp
 integer              :: ixam,iyam,izam
 integer              :: ixbm,iybm,izbm
 integer              :: mm
 real(dp),allocatable :: fmu(:)
 real(dp)             :: bigu
!=====

 ! Follow the notation of Obara and Saika, JCP  87 3963 (1986)
 zeta_ab = ga%alpha + gb%alpha
 ksi_ab   = ga%alpha * gb%alpha / zeta_ab
 ab2    = SUM( ( ga%x0(:)-gb%x0(:) )**2 )
 p(:)    = ( ga%alpha * ga%x0(:) + gb%alpha * gb%x0(:) ) / zeta_ab
 ap(:) =  p(:) - ga%x0(:)
 bp(:) =  p(:) - gb%x0(:)
 cp(:) =  p(:) - c(:)
 fact  = 0.5_dp / zeta_ab
 bigu  = zeta_ab * SUM( cp(:)**2 )


! v_tmp_x_m(:,:) =  100.0_dp
! v_tmp_y_m(:,:) =  100.0_dp
! v_tmp_z_m(:,:) =  100.0_dp

 v_tmp_x_mp1(:,:) =  0.0_dp
 v_tmp_y_mp1(:,:) =  0.0_dp
 v_tmp_z_mp1(:,:) =  0.0_dp

 do mm=ga%am+gb%am,0,-1
   allocate(fmu(0:mm))
   call boys_function_c(fmu,mm,bigu)

   !
   ! direction X
   !
   v_tmp_x_m(0,0) = 2.0 * fmu(mm) * (pi/zeta_ab) * EXP( - ksi_ab * ab2 )

   ixam=0
   ixbm=0
   do ix=1,ga%am+gb%am-mm
     do ixa=0,MIN(ga%nx,ix)
       ixb=ix-ixa
       if(ixb>gb%nx) cycle
       ixam=MAX(ixam,ixa)
       ixbm=MAX(ixbm,ixb)

       if(ixa>0) then
         ixap=ixa-1
         v_tmp_x_m(ixap+1,ixb) = ap(1) * v_tmp_x_m(ixap,ixb) - cp(1) * v_tmp_x_mp1(ixap,ixb)
         if(ixap>0) v_tmp_x_m(ixap+1,ixb) = v_tmp_x_m(ixap+1,ixb) + fact * ixap * ( v_tmp_x_m(ixap-1,ixb) -  v_tmp_x_mp1(ixap-1,ixb) )
         if(ixb>0)  v_tmp_x_m(ixap+1,ixb) = v_tmp_x_m(ixap+1,ixb) + fact * ixb  * ( v_tmp_x_m(ixap,ixb-1) -  v_tmp_x_mp1(ixap,ixb-1) )
       else
         ixbp=ixb-1
         v_tmp_x_m(ixa,ixbp+1) = bp(1) * v_tmp_x_m(ixa,ixbp) - cp(1) * v_tmp_x_mp1(ixa,ixbp)
         if(ixbp>0) v_tmp_x_m(ixa,ixbp+1) = v_tmp_x_m(ixa,ixbp+1) + fact * ixbp * ( v_tmp_x_m(ixa,ixbp-1) -  v_tmp_x_mp1(ixa,ixbp-1) )
         if(ixa>0)  v_tmp_x_m(ixa,ixbp+1) = v_tmp_x_m(ixa,ixbp+1) + fact * ixa  * ( v_tmp_x_m(ixa-1,ixbp) -  v_tmp_x_mp1(ixa-1,ixbp) )
       endif
  
     enddo
  
   enddo
  
   !
   ! direction Y
   !
   v_tmp_y_m(0,0) = v_tmp_x_m(ixam,ixbm)

   iyam=0
   iybm=0
   do iy=1,ga%am+gb%am-mm
     do iya=0,MIN(ga%ny,iy)
       iyb=iy-iya
       if(iyb>gb%ny) cycle
       iyam=MAX(iyam,iya)
       iybm=MAX(iybm,iyb)

       if(iya>0) then
         iyap=iya-1
         v_tmp_y_m(iyap+1,iyb) = ap(2) * v_tmp_y_m(iyap,iyb) - cp(2) * v_tmp_y_mp1(iyap,iyb)
         if(iyap>0) v_tmp_y_m(iyap+1,iyb) = v_tmp_y_m(iyap+1,iyb) + fact * iyap * ( v_tmp_y_m(iyap-1,iyb) -  v_tmp_y_mp1(iyap-1,iyb) )
         if(iyb>0)  v_tmp_y_m(iyap+1,iyb) = v_tmp_y_m(iyap+1,iyb) + fact * iyb  * ( v_tmp_y_m(iyap,iyb-1) -  v_tmp_y_mp1(iyap,iyb-1) )
       else
         iybp=iyb-1
         v_tmp_y_m(iya,iybp+1) = bp(2) * v_tmp_y_m(iya,iybp) - cp(2) * v_tmp_y_mp1(iya,iybp)
         if(iybp>0) v_tmp_y_m(iya,iybp+1) = v_tmp_y_m(iya,iybp+1) + fact * iybp * ( v_tmp_y_m(iya,iybp-1) -  v_tmp_y_mp1(iya,iybp-1) )
         if(iya>0)  v_tmp_y_m(iya,iybp+1) = v_tmp_y_m(iya,iybp+1) + fact * iya  * ( v_tmp_y_m(iya-1,iybp) -  v_tmp_y_mp1(iya-1,iybp) )
       endif
  
     enddo
  
   enddo
  
   !
   ! direction Z
   !
   v_tmp_z_m(0,0) = v_tmp_y_m(iyam,iybm)

   izam=0
   izbm=0
   do iz=1,ga%am+gb%am-mm
     do iza=0,MIN(ga%nz,iz)
       izb=iz-iza
       if(izb>gb%nz) cycle
       izam=MAX(izam,iza)
       izbm=MAX(izbm,izb)

       if(iza>0) then
         izap=iza-1
         v_tmp_z_m(izap+1,izb) = ap(3) * v_tmp_z_m(izap,izb) - cp(3) * v_tmp_z_mp1(izap,izb)
         if(izap>0) v_tmp_z_m(izap+1,izb) = v_tmp_z_m(izap+1,izb) + fact * izap * ( v_tmp_z_m(izap-1,izb) -  v_tmp_z_mp1(izap-1,izb) )
         if(izb>0)  v_tmp_z_m(izap+1,izb) = v_tmp_z_m(izap+1,izb) + fact * izb  * ( v_tmp_z_m(izap,izb-1) -  v_tmp_z_mp1(izap,izb-1) )
       else
         izbp=izb-1
         v_tmp_z_m(iza,izbp+1) = bp(3) * v_tmp_z_m(iza,izbp) - cp(3) * v_tmp_z_mp1(iza,izbp)
         if(izbp>0) v_tmp_z_m(iza,izbp+1) = v_tmp_z_m(iza,izbp+1) + fact * izbp * ( v_tmp_z_m(iza,izbp-1) -  v_tmp_z_mp1(iza,izbp-1) )
         if(iza>0)  v_tmp_z_m(iza,izbp+1) = v_tmp_z_m(iza,izbp+1) + fact * iza  * ( v_tmp_z_m(iza-1,izbp) -  v_tmp_z_mp1(iza-1,izbp) )
       endif
  
     enddo
  
   enddo

   v_tmp_x_mp1(0:ixam,0:ixbm) =  v_tmp_x_m(0:ixam,0:ixbm)
   v_tmp_y_mp1(0:iyam,0:iybm) =  v_tmp_y_m(0:iyam,0:iybm)
   v_tmp_z_mp1(0:izam,0:izbm) =  v_tmp_z_m(0:izam,0:izbm)

   deallocate(fmu)
 enddo


 v_ab = - zatom * v_tmp_z_m(izam,izbm) * ga%norm_factor * gb%norm_factor


end subroutine nucleus_recurrence