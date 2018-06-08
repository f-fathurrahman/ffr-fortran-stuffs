!=========================================================================
subroutine setup_nucleus_grad(print_matrix_,basis,hamiltonian_nucleus_grad)
  use m_definitions
 use m_atoms
  use m_timing
  use m_mpi
  use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: hamiltonian_nucleus_grad(basis%nbf,basis%nbf,natom+1,3)
!=====
 integer              :: ishell,jshell
 integer              :: ibf1,ibf2,jbf1,jbf2
 integer              :: natom_local
 integer              :: i_cart,j_cart
 integer              :: ni,nj,ni_cart,nj_cart,li,lj
 integer              :: iatom
 character(len=100)   :: title
 real(dp)             :: vnucleus_ij
 real(dp),allocatable :: matrixA(:,:)
 real(dp),allocatable :: matrixB(:,:)

 real(C_DOUBLE),allocatable        :: array_cart_gradAx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradAz(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBx(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBy(:)
 real(C_DOUBLE),allocatable        :: array_cart_gradBz(:)
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
 write(stdout,'(/,a)') ' Setup nucleus-electron part of the Hamiltonian gradient (LIBINT)'
 if( nproc_world > 1 ) then
   natom_local=0
   do iatom=1,natom
     if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle
     natom_local = natom_local + 1
   enddo
   write(stdout,'(a)')         '   Parallelizing over atoms'
   write(stdout,'(a,i5,a,i5)') '   this proc treats ',natom_local,' over ',natom
 endif

 hamiltonian_nucleus_grad(:,:,:,:) = 0.0_dp

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

     allocate(array_cart_gradAx(ni_cart*nj_cart))
     allocate(array_cart_gradAy(ni_cart*nj_cart))
     allocate(array_cart_gradAz(ni_cart*nj_cart))
     allocate(array_cart_gradBx(ni_cart*nj_cart))
     allocate(array_cart_gradBy(ni_cart*nj_cart))
     allocate(array_cart_gradBz(ni_cart*nj_cart))


     do iatom=1,natom
       if( rank_world /= MODULO(iatom-1,nproc_world) ) cycle


       C(:) = x(:,iatom)

#ifdef HAVE_LIBINT_ONEBODY
       call libint_elecpot_grad(amA,contrdepthA,A,alphaA,cA, &
                                amB,contrdepthB,B,alphaB,cB, &
                                C,                           &
                                array_cart_gradAx,array_cart_gradAy,array_cart_gradAz, &
                                array_cart_gradBx,array_cart_gradBy,array_cart_gradBz)
#else
       call die('nuclear potential gradient not implemented without LIBINT one-body terms')
#endif
       array_cart_gradAx(:) = array_cart_gradAx(:) * (-zatom(iatom))
       array_cart_gradAy(:) = array_cart_gradAy(:) * (-zatom(iatom))
       array_cart_gradAz(:) = array_cart_gradAz(:) * (-zatom(iatom))
       array_cart_gradBx(:) = array_cart_gradBx(:) * (-zatom(iatom))
       array_cart_gradBy(:) = array_cart_gradBy(:) * (-zatom(iatom))
       array_cart_gradBz(:) = array_cart_gradBz(:) * (-zatom(iatom))

       ! X
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAx,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBx,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,iatom  ,1) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,1) = hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,1) + matrixA(:,:)

       ! Y
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAy,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBy,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,iatom  ,2) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,2) = hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,2) + matrixA(:,:)

       ! Z
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradAz,matrixA)
       call transform_libint_to_molgw(basis%gaussian_type,li,lj,array_cart_gradBz,matrixB)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,iatom  ,3) = -matrixA(:,:) - matrixB(:,:)
       hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,3) = hamiltonian_nucleus_grad(ibf1:ibf2,jbf1:jbf2,natom+1,3) + matrixA(:,:)

     enddo
     deallocate(alphaA,cA)


     deallocate(array_cart_gradAx)
     deallocate(array_cart_gradAy)
     deallocate(array_cart_gradAz)
     deallocate(array_cart_gradBx)
     deallocate(array_cart_gradBy)
     deallocate(array_cart_gradBz)
     deallocate(matrixA)
     deallocate(matrixB)

   enddo
   deallocate(alphaB,cB)

 enddo

 !
 ! Reduce operation
 call xsum_world(hamiltonian_nucleus_grad)

 title='===  Nucleus potential contribution (LIBINT) C1X==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,1))
 title='===  Nucleus potential contribution (LIBINT) C1Y==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,2))
 title='===  Nucleus potential contribution (LIBINT) C1Z==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus_grad(:,:,1,3))

 call stop_clock(timing_hamiltonian_nuc)

end subroutine setup_nucleus_grad
