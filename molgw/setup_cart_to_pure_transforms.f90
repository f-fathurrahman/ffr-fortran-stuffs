!=========================================================================
subroutine setup_cart_to_pure_transforms()
  use m_definitions, only: dp, MOLGW_LMAX, stdout
  use m_cart_to_pure
 implicit none

!=====
 integer  :: ni,nic
 integer  :: ii,jj,kk
 integer  :: nx,ny,nz
 integer  :: il,im
 integer  :: it,iu,is
 real(dp) :: rtmp
!=====
  real(dp) :: double_factorial
  real(dp) :: cnk, ank
  integer :: number_basis_function_am

 write(stdout,'(/,1x,a,i2)') 'Setting up the cartesian to pure transforms up to l= ',MOLGW_LMAX

 allocate(cart_to_pure     (0:MOLGW_LMAX,2))
 allocate(cart_to_pure_norm(0:MOLGW_LMAX,2))

 !
 ! First setup trivial transforms in the case of CARTESIAN gaussians
 !
 do il=0,MOLGW_LMAX
   ni = number_basis_function_am('CART',il)
   allocate(cart_to_pure     (il,CARTG)%matrix(ni,ni))
   allocate(cart_to_pure_norm(il,CARTG)%matrix(ni,ni))
   cart_to_pure(il,CARTG)%matrix(:,:) = 0.0_dp
   do ii=1,ni
     cart_to_pure(il,CARTG)%matrix(ii,ii) = 1.0_dp
   enddo
 enddo


 !
 ! Second setup the complicated transforms in the case of PURE gaussians
 !
 do il=0,MOLGW_LMAX
   nic = number_basis_function_am('CART',il)
   ni  = number_basis_function_am('PURE',il)
   allocate(cart_to_pure(il,PUREG)%matrix(nic,ni))
   allocate(cart_to_pure_norm(il,PUREG)%matrix(nic,ni))
   cart_to_pure_norm(il,PUREG)%matrix(:,:) = 0.0_dp

   kk=0
   do ii=0,il
     nx = il - ii
     do jj=0,ii
       kk = kk + 1
       ny = ii - jj
       nz = jj

       rtmp = 0.0_dp
       do it=0,il/2
         do iu=0,it
           if( 2*it-2*iu == nx .AND. 2*iu == ny .AND. il - 2*it == nz ) then
             rtmp = rtmp + (-1)**it * 0.50_dp**(2*it) * cnk(il-it,it) * cnk(il,it) * cnk(it,iu)
           endif
         enddo
       enddo
       cart_to_pure_norm(il,PUREG)%matrix(kk,il+1) = rtmp / SQRT( double_factorial(2*il-1) )

       do im=1,il
         rtmp = 0.0_dp
         do it=0,(il-im)/2
           do iu=0,it
             do is=0,(im-1)/2
               if( im + 2*it - 2*iu - 2*is -1 == nx &
                 .AND. 2*iu + 2*is + 1 == ny &
                 .AND. il - im - 2*it == nz ) then
                 rtmp = rtmp + (-1)**(it+is) * 0.50_dp**(im+2*it) * SQRT( 2.0_dp * ank(il+im,il) / ank(il,il-im) ) &
                           * cnk(it,iu) * cnk(im,2*is+1) * cnk(il-it,im+it) * cnk(il,it)
               endif

             enddo
           enddo
         enddo
         cart_to_pure_norm(il,PUREG)%matrix(kk,il+1-im) = rtmp / SQRT( double_factorial(2*il-1) ) 
       enddo

       do im=1,il
         rtmp = 0.0_dp
         do it=0,(il-im)/2
           do iu=0,it
             do is=0,im/2
               if( im + 2*it - 2*iu - 2*is == nx &
                 .AND. 2*iu + 2*is == ny &
                 .AND. il - im - 2*it == nz ) then
                 rtmp = rtmp + (-1)**(it+is) * 0.50_dp**(im+2*it) * SQRT( 2.0_dp * ank(il+im,il) / ank(il,il-im) ) &
                           * cnk(it,iu) * cnk(im,2*is) * cnk(il-it,im+it) * cnk(il,it)
               endif

             enddo
           enddo
         enddo
         cart_to_pure_norm(il,PUREG)%matrix(kk,il+1+im) = rtmp / SQRT( double_factorial(2*il-1) )
       enddo


     enddo
   enddo
 enddo
 

 !
 ! Introduce the normalization coefficient part that depends on (nx,ny,nz)
 do il=0,MOLGW_LMAX
   kk=0
   do ii=0,il
     nx = il - ii
     do jj=0,ii
       kk = kk + 1
       ny = ii - jj
       nz = jj
       cart_to_pure_norm(il,CARTG)%matrix(kk,:) = cart_to_pure(il,CARTG)%matrix(kk,:) &
                 / SQRT( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) )

       cart_to_pure(il,PUREG)%matrix(kk,:) = cart_to_pure_norm(il,PUREG)%matrix(kk,:) &
                     * SQRT( double_factorial(2*nx-1) * double_factorial(2*ny-1) * double_factorial(2*nz-1) )

     enddo
   enddo
 enddo


 write(stdout,*) 'Transformations set up completed for both CARTESIAN and PURE Gaussians'
 write(stdout,*) 

end subroutine setup_cart_to_pure_transforms