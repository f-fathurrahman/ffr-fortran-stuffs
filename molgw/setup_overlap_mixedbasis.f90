!=========================================================================
subroutine setup_overlap_mixedbasis(print_matrix_,basis1,basis2,s_matrix)
  use m_definitions
  use m_basis_set
  use m_timing
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis1,basis2
 real(dp),intent(out)       :: s_matrix(basis1%nbf,basis2%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: ibf_cart,jbf_cart
 integer              :: i_cart,j_cart,ij
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: overlap

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

 call start_clock(timing_overlap)
 write(stdout,'(/,a)') ' Setup mixed overlap matrix S (LIBINT)'

 if( basis1%gaussian_type /= basis2%gaussian_type ) call die('setup_overlap_mixedbasis_libint: case not implemented')

 do jshell=1,basis2%nshell
   lj      = basis2%shell(jshell)%am
   nj_cart = number_basis_function_am('CART',lj)
   nj      = number_basis_function_am(basis2%gaussian_type,lj)
   jbf1    = basis2%shell(jshell)%istart
   jbf2    = basis2%shell(jshell)%iend

   call set_libint_shell(basis2%shell(jshell),amB,contrdepthB,B,alphaB,cB)

   do ishell=1,basis1%nshell
     li      = basis1%shell(ishell)%am
     ni_cart = number_basis_function_am('CART',li)
     ni      = number_basis_function_am(basis1%gaussian_type,li)
     ibf1    = basis1%shell(ishell)%istart
     ibf2    = basis1%shell(ishell)%iend

     call set_libint_shell(basis1%shell(ishell),amA,contrdepthA,A,alphaA,cA)


     allocate(array_cart(ni_cart*nj_cart))

#ifdef HAVE_LIBINT_ONEBODY

     call libint_overlap(amA,contrdepthA,A,alphaA,cA, &
                         amB,contrdepthB,B,alphaB,cB, &
                         array_cart)

     call transform_libint_to_molgw(basis1%gaussian_type,li,lj,array_cart,matrix)
#else
     ij = 0
     do i_cart=1,ni_cart
       do j_cart=1,nj_cart
         ij = ij + 1
         ibf_cart = basis1%shell(ishell)%istart_cart + i_cart - 1
         jbf_cart = basis2%shell(jshell)%istart_cart + j_cart - 1
         call overlap_basis_function(basis1%bfc(ibf_cart),basis2%bfc(jbf_cart),overlap)
         array_cart(ij) = overlap
       enddo
     enddo
     call transform_molgw_to_molgw(basis1%gaussian_type,li,lj,array_cart,matrix)
#endif

     deallocate(alphaA,cA)


     s_matrix(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)

     deallocate(array_cart,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 
 call stop_clock(timing_overlap)


end subroutine setup_overlap_mixedbasis
