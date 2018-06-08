!=========================================================================
subroutine setup_kinetic_grad(print_matrix_,basis,hamiltonian_kinetic_grad)
  use m_definitions
  use m_basis_set
  use m_timing
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic_grad(basis%nbf,basis%nbf,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradx(:)
 real(C_DOUBLE),allocatable        :: array_cart_grady(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradz(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
!=====
  integer :: number_basis_function_am

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup gradient of the kinetic part of the Hamiltonian (LIBINT)'

 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart_gradx(ni_cart*nj_cart))
     allocate(array_cart_grady(ni_cart*nj_cart))
     allocate(array_cart_gradz(ni_cart*nj_cart))

#ifdef HAVE_LIBINT_ONEBODY
     call libint_kinetic_grad(amA,contrdepthA,A,alphaA,cA, &
                              amB,contrdepthB,B,alphaB,cB, &
                              array_cart_gradx,array_cart_grady,array_cart_gradz)
#else
     call die('kinetic operator gradient not implemented without LIBINT one-body terms')
#endif
     deallocate(alphaA,cA)

     ! X
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradx,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,1) = matrix(:,:)

     ! Y
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_grady,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,2) = matrix(:,:)

     ! Z
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradz,matrix)
     hamiltonian_kinetic_grad(ibf1:ibf2,jbf1:jbf2,3) = matrix(:,:)


     deallocate(array_cart_gradx)
     deallocate(array_cart_grady)
     deallocate(array_cart_gradz)
     deallocate(matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='===  Kinetic energy contribution (LIBINT) X ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,1))
 title='===  Kinetic energy contribution (LIBINT) Y ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,2))
 title='===  Kinetic energy contribution (LIBINT) Z ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic_grad(:,:,3))

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic_grad
