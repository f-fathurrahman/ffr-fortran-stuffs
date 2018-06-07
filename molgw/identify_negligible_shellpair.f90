!=========================================================================
!
! Find negligible shell pairs with
! Cauchy-Schwarz inequality: (ij|1/r|kl)**2 <= (ij|1/r|ij) (kl|1/r|(kl) 
!
!=========================================================================
subroutine identify_negligible_shellpair(basis)
  USE m_definitions
  USE m_basis_set
  USE m_mpi
  USE m_eri
  USE m_timing
 implicit none

 type(basis_set),intent(in)   :: basis
!=====
 integer                      :: info,ip
 integer                      :: ibf,jbf
 integer                      :: n1c,n2c
 integer                      :: ni,nj
 integer                      :: ami,amj
 integer                      :: ishell,jshell
 real(dp),allocatable         :: integrals(:,:,:,:)
 real(dp)                     :: workload(nproc_world)
 integer                      :: shell_proc(basis%nshell)
!=====
! variables used to call C
 integer(C_INT)               :: am1,am2
 integer(C_INT)               :: ng1,ng2
 real(C_DOUBLE),allocatable   :: alpha1(:),alpha2(:)
 real(C_DOUBLE)               :: x01(3),x02(3)
 real(C_DOUBLE),allocatable   :: coeff1(:),coeff2(:)
 real(C_DOUBLE),allocatable   :: int_shell(:)
!=====
  INTEGER :: number_basis_function_am
  REAL(dp) :: cost_function_eri

 call start_clock(timing_eri_screening)
 write(stdout,'(/,a)')    ' Cauchy-Schwartz screening of the 3- or 4-center integrals'

 !
 ! Load balancing
 workload(:) = 0.0_dp
 do jshell=1,basis%nshell
   amj = basis%shell(jshell)%am
   ip = MINLOC(workload(:),DIM=1)
   !
   ! Cost function was evaluated from a few runs
   workload(ip) = workload(ip) + cost_function_eri(amj)
   shell_proc(jshell) = ip - 1
 enddo


 negligible_shellpair(:,:) = .TRUE.

 do jshell=1,basis%nshell
   !
   ! Workload is distributed here
   if( shell_proc(jshell) /= rank_world ) cycle

   amj = basis%shell(jshell)%am
   nj  = number_basis_function_am( basis%gaussian_type , amj )
   n2c = number_basis_function_am( 'CART' , amj )
   am2 = basis%shell(jshell)%am
   ng2 = basis%shell(jshell)%ng

   do ishell=1,basis%nshell
     ami = basis%shell(ishell)%am
     if( ami < amj ) cycle

     ni = number_basis_function_am( basis%gaussian_type , ami )
     n1c = number_basis_function_am( 'CART' , ami )
     am1 = basis%shell(ishell)%am
     ng1 = basis%shell(ishell)%ng

     allocate(alpha1(ng1),alpha2(ng2))
     allocate(coeff1(ng1),coeff2(ng2))
     alpha1(:) = basis%shell(ishell)%alpha(:)
     alpha2(:) = basis%shell(jshell)%alpha(:)
     x01(:) = basis%shell(ishell)%x0(:)
     x02(:) = basis%shell(jshell)%x0(:)
     coeff1(:) = basis%shell(ishell)%coeff(:)
     coeff2(:) = basis%shell(jshell)%coeff(:)

     allocate( int_shell( n1c*n2c*n1c*n2c ) )


     call libint_4center(am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         am1,ng1,x01,alpha1,coeff1, &
                         am2,ng2,x02,alpha2,coeff2, &
                         0.0_C_DOUBLE,int_shell)

     call transform_libint_to_molgw(basis%gaussian_type,ami,amj,ami,amj,int_shell,integrals)


     do ibf=1,ni
       do jbf=1,nj
         if( ABS( integrals(ibf,jbf,ibf,jbf) ) > TOL_INT**2 ) negligible_shellpair(ishell,jshell) = .FALSE.
       enddo
     enddo

     !
     ! Symmetrize
     negligible_shellpair(jshell,ishell) = negligible_shellpair(ishell,jshell)

     deallocate(integrals)
     deallocate(int_shell)
     deallocate(alpha1,alpha2)
     deallocate(coeff1,coeff2)

   enddo
 enddo

 call xand_world(negligible_shellpair)

 call stop_clock(timing_eri_screening)


end subroutine identify_negligible_shellpair
