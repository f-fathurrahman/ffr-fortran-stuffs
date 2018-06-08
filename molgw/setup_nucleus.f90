!=========================================================================
subroutine setup_nucleus(print_matrix_,basis,hamiltonian_nucleus)
 use m_definitions
 use m_atoms
 use m_basis_set
 use m_mpi
 use m_timing
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: natom_local
 integer              :: ibf_cart,jbf_cart,ij
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp),allocatable :: matrix(:,:)
 real(dp)             :: nucleus

 real(C_DOUBLE),allocatable        :: array_cart(:)
 real(C_DOUBLE),allocatable        :: array_cart_C(:)
 integer(C_INT)                    :: amA,contrdepthA
 real(C_DOUBLE)                    :: A(3)
 real(C_DOUBLE),allocatable        :: alphaA(:)
 real(C_DOUBLE),allocatable        :: cA(:)
 integer(C_INT)                    :: amB,contrdepthB
 real(C_DOUBLE)                    :: B(3)
 real(C_DOUBLE),allocatable        :: alphaB(:)
 real(C_DOUBLE),allocatable        :: cB(:)
 real(C_DOUBLE)                    :: C(3)
!=====
  integer :: number_basis_function_am

 call start_clock(timing_hamiltonian_nuc)
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian (LIBINT)'
 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif


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
     allocate(array_cart_C(ni_cart*nj_cart))
     array_cart(:) = 0.0_dp

     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle

       C(:) = x(:,iatom)
#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot(amA,contrdepthA,A,alphaA,cA, &
                           amB,contrdepthB,B,alphaB,cB, &
                           C,array_cart_C)
       array_cart(:) = array_cart(:) - zatom(iatom) * array_cart_C(:) 
#else
       ij = 0
       do i_cart=1,ni_cart
         do j_cart=1,nj_cart
           ij = ij + 1
           ibf_cart = basis%shell(ishell)%istart_cart + i_cart - 1
           jbf_cart = basis%shell(jshell)%istart_cart + j_cart - 1
           call nucleus_basis_function(basis%bfc(ibf_cart),basis%bfc(jbf_cart),zatom(iatom),x(:,iatom),nucleus)
           array_cart(ij) = array_cart(ij) + nucleus
         enddo
       enddo
#endif

     enddo
     deallocate(alphaA,cA)

#ifdef HAVE_LIBINT_ONEBODY
     call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#else
     call transform_molgw_to_molgw(basis%gaussian_type,li,lj,array_cart,matrix)
#endif

     hamiltonian_nucleus(ibf1:ibf2,jbf1:jbf2) = matrix(:,:)
     hamiltonian_nucleus(jbf1:jbf2,ibf1:ibf2) = TRANSPOSE(matrix(:,:))


     deallocate(array_cart,array_cart_C,matrix)

   enddo
   deallocate(alphaB,cB)
 enddo

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus)

 title='===  Nucleus potential contribution (LIBINT) ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus
