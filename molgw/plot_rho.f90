!=========================================================================
subroutine plot_rho(nstate,basis,occupation,c_matrix)
 use m_definitions
 use m_mpi
 use m_atoms
 use m_cart_to_pure
 use m_inputparam, only: nspin
 use m_basis_set
 use m_dft_grid,only: calculate_basis_functions_r
 implicit none
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer,parameter          :: nr=5000
 real(dp),parameter         :: length=4.0_dp
 integer                    :: gt
 integer                    :: ir,ibf
 integer                    :: istate1,istate2,istate,ispin
 real(dp)                   :: rr(3)
 real(dp),allocatable       :: phi(:,:)
 real(dp)                   :: u(3),a(3)
 logical                    :: file_exists
 real(dp)                   :: xxmin,xxmax
 real(dp)                   :: basis_function_r(basis%nbf)
 integer                    :: rhorfile
!=====

 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)') 'Plotting the density'

 gt = get_gaussian_type_tag(basis%gaussian_type)

 inquire(file='manual_plotrho',exist=file_exists)
 if(file_exists) then
   open(newunit=rhorfile,file='manual_plotrho',status='old')
   read(rhorfile,*) u(:)
   read(rhorfile,*) a(:)
   close(rhorfile)
 else
   u(:)=0.0_dp
   u(1)=1.0_dp
   a(:)=0.0_dp
 endif
 u(:) = u(:) / NORM2(u)
 allocate(phi(nstate,nspin))
 write(stdout,'(a,3(2x,f8.3))') ' direction:',u(:)
 write(stdout,'(a,3(2x,f8.3))') ' origin:   ',a(:)

 xxmin = MINVAL( u(1)*x(1,:) + u(2)*x(2,:) + u(3)*x(3,:) ) - length
 xxmax = MAXVAL( u(1)*x(1,:) + u(2)*x(2,:) + u(3)*x(3,:) ) + length


 do ir=1,nr
   rr(:) = ( xxmin + (ir-1)*(xxmax-xxmin)/REAL(nr-1,dp) ) * u(:) + a(:)


   call calculate_basis_functions_r(basis,rr,basis_function_r)

   do ispin=1,nspin
     phi(:,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,:,ispin) )
   enddo

   write(103,'(50(e16.8,2x))') DOT_PRODUCT(rr(:),u(:)),SUM( phi(:,:)**2 * occupation(:,:) )

   write(104,'(50(e16.8,2x))') NORM2(rr(:)),SUM( phi(:,:)**2 * occupation(:,:) ) * 4.0_dp * pi * NORM2(rr(:))**2
!   write(105,'(50(e16.8,2x))') NORM2(rr(:)),SUM( phi(1,:)**2 * occupation(1,:) ) * 4.0_dp * pi * NORM2(rr(:))**2
!   write(106,'(50(e16.8,2x))') NORM2(rr(:)),SUM( phi(2:,:)**2 * occupation(2:,:) ) * 4.0_dp * pi * NORM2(rr(:))**2

 enddo

 deallocate(phi)

end subroutine plot_rho