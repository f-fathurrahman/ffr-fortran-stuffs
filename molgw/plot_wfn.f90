!=========================================================================
subroutine plot_wfn(nstate,basis,c_matrix)
 use m_definitions
 use m_mpi
 use m_inputparam, only: nspin
 use m_atoms
 use m_cart_to_pure
 use m_basis_set
 use m_dft_grid,only: calculate_basis_functions_r
 implicit none
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer,parameter          :: nr=2000
 real(dp),parameter         :: length=10.0_dp
 integer                    :: gt
 integer                    :: ir,ibf
 integer                    :: istate1,istate2,istate,ispin
 real(dp)                   :: rr(3)
 real(dp),allocatable       :: phi(:,:),phase(:,:)
 real(dp)                   :: u(3),a(3)
 logical                    :: file_exists
 real(dp)                   :: xxmin,xxmax
 real(dp)                   :: basis_function_r(basis%nbf)
 integer                    :: wfrfile
!=====

 if( .NOT. is_iomaster ) return

 gt = get_gaussian_type_tag(basis%gaussian_type)

 write(stdout,'(/,1x,a)') 'Plotting some selected wavefunctions'
 inquire(file='manual_plotwfn',exist=file_exists)
 if(file_exists) then
   open(newunit=wfrfile,file='manual_plotwfn',status='old')
   read(wfrfile,*) istate1,istate2
   read(wfrfile,*) u(:)
   read(wfrfile,*) a(:)
   close(wfrfile)
 else
   istate1=1
   istate2=2
   u(:)=0.0_dp
   u(1)=1.0_dp
   a(:)=0.0_dp
 endif
 u(:) = u(:) / SQRT(SUM(u(:)**2))
 allocate(phase(istate1:istate2,nspin),phi(istate1:istate2,nspin))
 write(stdout,'(a,2(2x,i4))')   ' states:   ',istate1,istate2
 write(stdout,'(a,3(2x,f8.3))') ' direction:',u(:)
 write(stdout,'(a,3(2x,f8.3))') ' origin:   ',a(:)

 xxmin = MINVAL( u(1)*x(1,:) + u(2)*x(2,:) + u(3)*x(3,:) ) - length
 xxmax = MAXVAL( u(1)*x(1,:) + u(2)*x(2,:) + u(3)*x(3,:) ) + length

 phase(:,:)=1.0_dp

 do ir=1,nr
   rr(:) = ( xxmin + (ir-1)*(xxmax-xxmin)/REAL(nr-1,dp) ) * u(:) + a(:)

   
   call calculate_basis_functions_r(basis,rr,basis_function_r)

   do ispin=1,nspin
     phi(istate1:istate2,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,istate1:istate2,ispin) )
   enddo

   !
   ! turn the wfns so that they are all positive at a given point
   if(ir==1) then
     do ispin=1,nspin
       do istate=istate1,istate2
         if( phi(istate,ispin) < 0.0_dp ) phase(istate,ispin) = -1.0_dp
       enddo
     enddo
   endif

   write(101,'(50(e16.8,2x))') DOT_PRODUCT(rr(:),u(:)),phi(:,:)*phase(:,:)
   write(102,'(50(e16.8,2x))') DOT_PRODUCT(rr(:),u(:)),phi(:,:)**2

 enddo

 deallocate(phase,phi)

end subroutine plot_wfn