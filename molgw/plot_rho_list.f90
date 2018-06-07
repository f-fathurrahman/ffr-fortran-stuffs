!=========================================================================
subroutine plot_rho_list(nstate,basis,occupation,c_matrix)
 use m_definitions
 use m_mpi
 use m_atoms
 use m_cart_to_pure
 use m_inputparam, only: nspin
 use m_basis_set
 implicit none
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
!=====
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
 integer                    :: ix,iy,iz
 integer,parameter          :: nx=75 ! 87
 integer,parameter          :: ny=75 ! 91
 integer,parameter          :: nz=90 ! 65
 real(dp),parameter         :: dx = 0.174913 ! 0.204034
 real(dp)                   :: rr0(3)
 integer                    :: unitfile
!=====
  INTEGER :: get_gaussian_type_tag

 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)') 'Plotting the density'

 gt = get_gaussian_type_tag(basis%gaussian_type)

 inquire(file='manual_plotrho',exist=file_exists)
 if(file_exists) then
   open(newunit=rhorfile,file='manual_plotrho',status='old')
   close(rhorfile)
 else
 endif
 allocate(phi(nstate,nspin))

 rr0(1) = -6.512752 ! -8.790885
 rr0(2) = -6.512752 ! -9.143313 
 rr0(3) = -7.775444 ! -6.512752

 open(newunit=unitfile,file='rho.dat',action='WRITE')
 do ix=1,nx
 do iy=1,ny
 do iz=1,nz
   rr(1) = ix-1
   rr(2) = iy-1
   rr(3) = iz-1
   rr(:) = rr0(:) + rr(:) * dx 


   call calculate_basis_functions_r(basis,rr,basis_function_r)

   do ispin=1,nspin
     phi(:,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,:,ispin) )
   enddo

   write(unitfile,'(1x,e16.8)') SUM( phi(3:,:)**2 * occupation(3:,:) )

 enddo
 enddo
 enddo
 close(unitfile)

 deallocate(phi)

end subroutine plot_rho_list
