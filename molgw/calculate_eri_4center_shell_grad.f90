!=========================================================================
subroutine calculate_eri_4center_shell_grad(basis,rcut,ijshellpair,klshellpair,&
                                            shell_gradA,shell_gradB,shell_gradC,shell_gradD)
  use m_definitions
  use m_basis_set
  use m_eri
  implicit none
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
 integer,intent(in)           :: ijshellpair,klshellpair
 real(dp),intent(out)         :: shell_gradA(:,:,:,:,:)
 real(dp),intent(out)         :: shell_gradB(:,:,:,:,:)
 real(dp),intent(out)         :: shell_gradC(:,:,:,:,:)
 real(dp),intent(out)         :: shell_gradD(:,:,:,:,:)
!=====
 logical                      :: is_longrange
 integer                      :: ishell,jshell,kshell,lshell
 integer                      :: n1c,n2c,n3c,n4c
 integer                      :: ni,nj,nk,nl
 integer                      :: ami,amj,amk,aml
 integer                      :: ibf,jbf,kbf,lbf
 real(dp),allocatable         :: grad_tmp(:,:,:,:)
!=====
! variables used to call C
 real(C_DOUBLE)               :: rcut_libint
 integer(C_INT)               :: ng1,ng2,ng3,ng4
 integer(C_INT)               :: am1,am2,am3,am4
 real(C_DOUBLE)               :: x01(3),x02(3),x03(3),x04(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:),coeff3(:),coeff4(:)
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:),alpha3(:),alpha4(:)
 real(C_DOUBLE),allocatable   :: gradAx(:)
 real(C_DOUBLE),allocatable   :: gradAy(:)
 real(C_DOUBLE),allocatable   :: gradAz(:)
 real(C_DOUBLE),allocatable   :: gradBx(:)
 real(C_DOUBLE),allocatable   :: gradBy(:)
 real(C_DOUBLE),allocatable   :: gradBz(:)
 real(C_DOUBLE),allocatable   :: gradCx(:)
 real(C_DOUBLE),allocatable   :: gradCy(:)
 real(C_DOUBLE),allocatable   :: gradCz(:)
 real(C_DOUBLE),allocatable   :: gradDx(:)
 real(C_DOUBLE),allocatable   :: gradDy(:)
 real(C_DOUBLE),allocatable   :: gradDz(:)
!=====
  integer :: number_basis_function_am

 is_longrange = (rcut > 1.0e-12_dp)
 rcut_libint = rcut

 kshell = index_shellpair(1,klshellpair)
 lshell = index_shellpair(2,klshellpair)

 !
 ! The angular momenta are already ordered so that libint is pleased
 ! 1) amk+aml >= ami+amj
 ! 2) amk>=aml
 ! 3) ami>=amj
 amk = basis%shell(kshell)%am
 aml = basis%shell(lshell)%am


 ishell = index_shellpair(1,ijshellpair)
 jshell = index_shellpair(2,ijshellpair)

 ami = basis%shell(ishell)%am
 amj = basis%shell(jshell)%am
 if( amk+aml < ami+amj ) call die('calculate_4center_shell_grad: wrong ordering')

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

 allocate(gradAx(n1c*n2c*n3c*n4c))
 allocate(gradAy(n1c*n2c*n3c*n4c))
 allocate(gradAz(n1c*n2c*n3c*n4c))
 allocate(gradBx(n1c*n2c*n3c*n4c))
 allocate(gradBy(n1c*n2c*n3c*n4c))
 allocate(gradBz(n1c*n2c*n3c*n4c))
 allocate(gradCx(n1c*n2c*n3c*n4c))
 allocate(gradCy(n1c*n2c*n3c*n4c))
 allocate(gradCz(n1c*n2c*n3c*n4c))
 allocate(gradDx(n1c*n2c*n3c*n4c))
 allocate(gradDy(n1c*n2c*n3c*n4c))
 allocate(gradDz(n1c*n2c*n3c*n4c))

#ifdef HAVE_LIBINT_ONEBODY
 call libint_4center_grad(am1,ng1,x01,alpha1,coeff1, &
                          am2,ng2,x02,alpha2,coeff2, &
                          am3,ng3,x03,alpha3,coeff3, &
                          am4,ng4,x04,alpha4,coeff4, &
                          rcut_libint, &
                          gradAx,gradAy,gradAz, &
                          gradBx,gradBy,gradBz, &
                          gradCx,gradCy,gradCz, &
                          gradDx,gradDy,gradDz)
#endif

 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradAx,grad_tmp)
 shell_gradA(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradAy,grad_tmp)
 shell_gradA(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradAz,grad_tmp)
 shell_gradA(:,:,:,:,3) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradBx,grad_tmp)
 shell_gradB(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradBy,grad_tmp)
 shell_gradB(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradBz,grad_tmp)
 shell_gradB(:,:,:,:,3) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradCx,grad_tmp)
 shell_gradC(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradCy,grad_tmp)
 shell_gradC(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradCz,grad_tmp)
 shell_gradC(:,:,:,:,3) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradDx,grad_tmp)
 shell_gradD(:,:,:,:,1) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradDy,grad_tmp)
 shell_gradD(:,:,:,:,2) = grad_tmp(:,:,:,:) 
 call transform_libint_to_molgw(basis%gaussian_type,ami,amj,amk,aml,gradDz,grad_tmp)
 shell_gradD(:,:,:,:,3) = grad_tmp(:,:,:,:) 


 deallocate(grad_tmp)
 deallocate(gradAx)
 deallocate(gradAy)
 deallocate(gradAz)
 deallocate(gradBx)
 deallocate(gradBy)
 deallocate(gradBz)
 deallocate(gradCx)
 deallocate(gradCy)
 deallocate(gradCz)
 deallocate(gradDx)
 deallocate(gradDy)
 deallocate(gradDz)
 deallocate(alpha1,alpha2,alpha3,alpha4)
 deallocate(coeff1)
 deallocate(coeff2)
 deallocate(coeff3)
 deallocate(coeff4)



end subroutine calculate_eri_4center_shell_grad
