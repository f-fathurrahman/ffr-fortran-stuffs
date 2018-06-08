!=========================================================================
subroutine setup_kinetic(print_matrix_,basis,hamiltonian_kinetic)
  use m_definitions
  use m_basis_set
  use m_timing
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_kinetic(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,ij
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: kinetic

 real(C_DOUBLE),allocatable :: array_cart(:)
 integer(C_INT)             :: amA,contrdepthA
 real(C_DOUBLE)             :: A(3)
 real(C_DOUBLE),allocatable :: alphaA(:)
 real(C_DOUBLE),allocatable :: cA(:)
 integer(C_INT)             :: amB,contrdepthB
 real(C_DOUBLE)             :: B(3)
 real(C_DOUBLE),allocatable :: alphaB(:)
 real(C_DOUBLE),allocatable :: cB(:)
!=====
  integer :: number_basis_function_am

 call start_clock(timing_hamiltonian_kin)
 write(stdout,'(/,a)') ' Setup kinetic part of the Hamiltonian (LIBINT)'


 do jshell=1,basis%nshell
   lj      = basis%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis%gaussian_type,lj)
   jbf1    = basis%shell(jshell)%istart
   jbf2    = basis%shell(jshell)%iend

   call set_libint_shell(basis%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=jshell,basis%nshell
     li      = basis%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis%gaussian_type,li)
     ibf1    = basis%shell(ishell)%istart
     ibf2    = basis%shell(ishell)%iend

     call set_libint_shell(basis%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))


#ifdef HAVE_LIBINT_ONEBODY
     call libint_kinetic(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
         call kinetic_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),kinetic)
         array_cart(ij) = kinetic
       enddo
     enddo
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif
     deallocate(alphaA,cA)



     hamiltonian_kinetic(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     hamiltonian_kinetic(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))


     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 title='===  Kinetic energy contribution (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_kinetic)

 call stop_clock(timing_hamiltonian_kin)

end subroutine setup_kinetic
