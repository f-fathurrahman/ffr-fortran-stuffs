!=========================================================================
subroutine calculate_eri_4center(basis,rcut)
  use m_definitions
  use m_basis_set
  use m_eri
 implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!=====
 logical                      :: is_longrange
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: ijshellpair,klshellpair
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
 real(dp),allocatable         :: integrals(:,:,:,:)
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====
  integer :: number_basis_function_am
  integer :: index_eri

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut
 if( .NOT. is_longrange ) then 
   write(stdout,'(/,a)') ' Calculate and store the 4-center Coulomb integrals'
 else
   write(stdout,'(/,a)') ' Calculate and store the 4-center Long-Range-only Coulomb integrals'
 endif



 do klshellpair=1,nshellpair
   kshell = index_shellpair(1,klshellpair)
   lshell = index_shellpair(2,klshellpair)

   !
   ! The angular momenta are already ordered so that libint is pleased
   ! 1) amk+aml >= ami+amj
   ! 2) amk>=aml
   ! 3) ami>=amj
   amk = basis%shell(kshell)%am
   aml = basis%shell(lshell)%am


   do ijshellpair=1,nshellpair
     ishell = index_shellpair(1,ijshellpair)
     jshell = index_shellpair(2,ijshellpair)

     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     if( amk+aml < ami+amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     nj = number_basis_function_am( basis%gaussian_type , amj )
     nk = number_basis_function_am( basis%gaussian_type , amk )
     nl = number_basis_function_am( basis%gaussian_type , aml )


     am1 = basis%shell(ishell)%am
     am2 = basis%shell(jshell)%am
     am3 = basis%shell(kshell)%am
     am4 = basis%shell(lshell)%am
     n1c = number_basis_function_am( 'CART' , ami )
     n2c = number_basis_function_am( 'CART' , amj )
     n3c = number_basis_function_am( 'CART' , amk )
     n4c = number_basis_function_am( 'CART' , aml )
     ng1 = basis%shell(ishell)%ng
     ng2 = basis%shell(jshell)%ng
     ng3 = basis%shell(kshell)%ng
     ng4 = basis%shell(lshell)%ng
     allocate(alpha1(ng1),alpha2(ng2),alpha3(ng3),alpha4(ng4))
     alpha1(:) = basis%shell(ishell)%alpha(:) 
     alpha2(:) = basis%shell(jshell)%alpha(:)
     alpha3(:) = basis%shell(kshell)%alpha(:)
     alpha4(:) = basis%shell(lshell)%alpha(:)
     x01(:) = basis%shell(ishell)%x0(:)
     x02(:) = basis%shell(jshell)%x0(:)
     x03(:) = basis%shell(kshell)%x0(:)
     x04(:) = basis%shell(lshell)%x0(:)
     allocate(coeff1(basis%shell(ishell)%ng))
     allocate(coeff2(basis%shell(jshell)%ng))
     allocate(coeff3(basis%shell(kshell)%ng))
     allocate(coeff4(basis%shell(lshell)%ng))
     coeff1(:)=basis%shell(ishell)%coeff(:)
     coeff2(:)=basis%shell(jshell)%coeff(:)
     coeff3(:)=basis%shell(kshell)%coeff(:)
     coeff4(:)=basis%shell(lshell)%coeff(:)

     allocate(int_shell(n1c*n2c*n3c*n4c))

     call libint_4center(am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         am3,ng3,x03,alpha3,coeff3, &
                         am4,ng4,x04,alpha4,coeff4, &
                         rcut_libint,int_shell)

     call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,int_shell,integrals)

     if( .NOT. is_longrange ) then
       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               eri_4center( index_eri(basis%shell(ishell)%istart+ibf-1, &
                                      basis%shell(jshell)%istart+jbf-1, &
                                      basis%shell(kshell)%istart+kbf-1, &
                                      basis%shell(lshell)%istart+lbf-1) ) = integrals(ibf,jbf,kbf,lbf)
             enddo
           enddo
         enddo
       enddo
     else
       do lbf=1,nl
         do kbf=1,nk
           do jbf=1,nj
             do ibf=1,ni
               eri_4center_lr( index_eri(basis%shell(ishell)%istart+ibf-1, &
                                         basis%shell(jshell)%istart+jbf-1, &
                                         basis%shell(kshell)%istart+kbf-1, &
                                         basis%shell(lshell)%istart+lbf-1) ) = integrals(ibf,jbf,kbf,lbf)
             enddo
           enddo
         enddo
       enddo
     endif


     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,alpha2,alpha3,alpha4)
     deallocate(coeff1)
     deallocate(coeff2)
     deallocate(coeff3)
     deallocate(coeff4)


   enddo
 enddo


 write(stdout,'(a,/)') ' All ERI have been calculated'


end subroutine calculate_eri_4center
