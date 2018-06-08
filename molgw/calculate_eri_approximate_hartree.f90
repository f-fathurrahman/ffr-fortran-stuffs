!=========================================================================
subroutine calculate_eri_approximate_hartree(basis,mv,nv,x0_rho,ng_rho,coeff_rho,alpha_rho,vhrho)
  use m_definitions
  use m_basis_set
  use m_scalapack
  use m_eri
  use m_cart_to_pure
 implicit none
 type(basis_set),intent(in)   :: basis
 integer,intent(in)           :: mv,nv
 real(dp),intent(in)          :: x0_rho(3)
 integer,intent(in)           :: ng_rho
 real(dp),intent(in)          :: coeff_rho(ng_rho),alpha_rho(ng_rho)
 real(dp),intent(inout)       :: vhrho(mv,nv)
!=====
 integer                      :: kshell,lshell
 integer                      :: klshellpair
 integer                      :: n3c,n4c
 integer                      :: nk,nl
 integer                      :: amk,aml
 integer                      :: kbf,lbf
 real(dp),allocatable         :: integrals(:,:,:)
 integer                      :: ilocal,jlocal,iglobal,jglobal
!=====
! variables used to call C
 integer(C_INT)               :: am1,am3,am4
 integer(C_INT)               :: ng1,ng3,ng4
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha3(:),alpha4(:)
 real(C_DOUBLE)               :: x01(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====
  integer :: number_basis_function_am

 if( parallel_ham .AND. parallel_buffer ) then    
   if( mv /= basis%nbf .OR. nv /= basis%nbf ) call die('calculate_eri_approximate_hartree: wrong dimension for the buffer')
 endif

 ! Nullify vhrho just for safety.
 ! I guess this is useless.
 if( .NOT. ( parallel_ham .AND. parallel_buffer )  ) then    
   vhrho(:,:) = 0.0_dp
 endif

 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am

   nk = number_basis_function_am( basis%gaussian_type , amk )
   nl = number_basis_function_am( basis%gaussian_type , aml )

   am1 = 0
   am3 = amk
   am4 = aml
   n3c = number_basis_function_am( 'CART' , amk )
   n4c = number_basis_function_am( 'CART' , aml )
   ng1 = ng_rho
   ng3 = basis%shell(kshell)%ng
   ng4 = basis%shell(lshell)%ng
   allocate(alpha1(ng1),alpha3(ng3),alpha4(ng4))
   allocate(coeff1(ng1),coeff3(ng3),coeff4(ng4))
   alpha1(:) = alpha_rho(:)
   alpha3(:) = basis%shell(kshell)%alpha(:)
   alpha4(:) = basis%shell(lshell)%alpha(:)
   coeff1(:) = coeff_rho(:) / 2.0_dp**1.25_dp / pi**0.75_dp * alpha_rho(:)**1.5_dp * cart_to_pure_norm(0,PUREG)%matrix(1,1)
   coeff3(:) = basis%shell(kshell)%coeff(:)
   coeff4(:) = basis%shell(lshell)%coeff(:)
   x01(:) = x0_rho(:)
   x03(:) = basis%shell(kshell)%x0(:)
   x04(:) = basis%shell(lshell)%x0(:)

   allocate( int_shell(n3c*n4c) )

   call libint_3center(am1,ng1,x01,alpha1,coeff1, &
                       am3,ng3,x03,alpha3,coeff3, &
                       am4,ng4,x04,alpha4,coeff4, &
                       0.0_C_DOUBLE,int_shell)

   call transform_libint_to_molgw(basis%gaussian_type,0,basis%gaussian_type,amk,aml,int_shell,integrals)

     

   if( parallel_ham .AND. parallel_buffer ) then    
     do lbf=1,nl
       do kbf=1,nk
         iglobal = basis%shell(kshell)%istart+kbf-1
         jglobal = basis%shell(lshell)%istart+lbf-1
         ilocal = iglobal
         jlocal = jglobal
         if( kshell == lshell ) then ! To avoid double-counting   
           vhrho(ilocal,jlocal) = vhrho(ilocal,jlocal) + integrals(1,kbf,lbf)  * 0.5_dp 
         else
           vhrho(ilocal,jlocal) = vhrho(ilocal,jlocal) + integrals(1,kbf,lbf) 
         endif
       enddo
     enddo

   else

     do lbf=1,nl
       do kbf=1,nk
         iglobal = basis%shell(kshell)%istart+kbf-1
         jglobal = basis%shell(lshell)%istart+lbf-1
         ilocal = rowindex_global_to_local('H',iglobal)
         jlocal = colindex_global_to_local('H',jglobal)
         if( ilocal*jlocal /= 0 ) then
           vhrho(ilocal,jlocal) = integrals(1,kbf,lbf)
         endif
         ! And the symmetric too
         iglobal = basis%shell(lshell)%istart+lbf-1
         jglobal = basis%shell(kshell)%istart+kbf-1
         ilocal = rowindex_global_to_local('H',iglobal)
         jlocal = colindex_global_to_local('H',jglobal)
         if( ilocal*jlocal /= 0 ) then
           vhrho(ilocal,jlocal) = integrals(1,kbf,lbf)
         endif
       enddo
     enddo
   endif


   deallocate(integrals)
   deallocate(int_shell)
   deallocate(alpha1,alpha3,alpha4)
   deallocate(coeff1,coeff3,coeff4)

 enddo


end subroutine calculate_eri_approximate_hartree
