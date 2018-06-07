!=========================================================================
subroutine plot_cube_wfn(nstate,basis,occupation,c_matrix)
 use m_definitions
 use m_mpi
 use m_inputparam, only: nspin,spin_fact
 use m_cart_to_pure
 use m_atoms
 use m_basis_set
 use m_dft_grid,only: calculate_basis_functions_r
 implicit none
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: occupation(nstate,nspin)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
!=====
 integer                    :: gt
 integer                    :: nx,ny,nz
 real(dp),parameter         :: length=3.499470_dp
 integer                    :: istate1,istate2,istate,ispin
 real(dp)                   :: rr(3)
 real(dp),allocatable       :: phi(:,:)
 real(dp)                   :: u(3),a(3)
 logical                    :: file_exists
 real(dp)                   :: xmin,xmax,ymin,ymax,zmin,zmax
 real(dp)                   :: dx,dy,dz
 real(dp)                   :: basis_function_r(basis%nbf)
 integer                    :: ix,iy,iz,iatom
 integer,allocatable        :: ocubefile(:,:)
 integer                    :: ocuberho(nspin)
 character(len=200)         :: file_name
 integer                    :: icubefile
!=====

 if( .NOT. is_iomaster ) return

 write(stdout,'(/,1x,a)') 'Plotting some selected wavefunctions in a cube file'

 gt = get_gaussian_type_tag(basis%gaussian_type)

 inquire(file='manual_cubewfn',exist=file_exists)
 if(file_exists) then
   open(newunit=icubefile,file='manual_cubewfn',status='old')
   read(icubefile,*) istate1,istate2
   read(icubefile,*) nx,ny,nz
   close(icubefile)
 else
   istate1=1
   istate2=2
   nx=40
   ny=40
   nz=40
 endif
 allocate(phi(istate1:istate2,nspin))
 write(stdout,'(a,2(2x,i4))')   ' states:   ',istate1,istate2

 xmin = MINVAL( x(1,:) ) - length
 xmax = MAXVAL( x(1,:) ) + length
 ymin = MINVAL( x(2,:) ) - length
 ymax = MAXVAL( x(2,:) ) + length
 zmin = MINVAL( x(3,:) ) - length
 zmax = MAXVAL( x(3,:) ) + length
 dx = (xmax-xmin)/REAL(nx,dp)
 dy = (ymax-ymin)/REAL(ny,dp)
 dz = (zmax-zmin)/REAL(nz,dp)
! xmin = -15.001591d0
! ymin = -15.001591d0
! zmin = -17.037892d0
! dx = 0.262502_dp
! dy = 0.262502_dp
! dz = 0.262502_dp

 allocate(ocubefile(istate1:istate2,nspin))

 do istate=istate1,istate2
   do ispin=1,nspin
     write(file_name,'(a,i3.3,a,i1,a)') 'wfn_',istate,'_',ispin,'.cube'
     open(newunit=ocubefile(istate,ispin),file=file_name)
     write(ocubefile(istate,ispin),'(a)') 'cube file generated from MOLGW'
     write(ocubefile(istate,ispin),'(a,i4)') 'wavefunction ',istate1
     write(ocubefile(istate,ispin),'(i6,3(f12.6,2x))') natom,xmin,ymin,zmin
     write(ocubefile(istate,ispin),'(i6,3(f12.6,2x))') nx,dx,0.,0.
     write(ocubefile(istate,ispin),'(i6,3(f12.6,2x))') ny,0.,dy,0.
     write(ocubefile(istate,ispin),'(i6,3(f12.6,2x))') nz,0.,0.,dz
     do iatom=1,natom
       write(ocubefile(istate,ispin),'(i6,4(2x,f12.6))') basis_element(iatom),0.0,x(:,iatom)
     enddo
   enddo
 enddo

 !
 ! check whether istate1:istate2 spans all the occupied states
 if( ALL( occupation(istate2+1,:) < completely_empty ) ) then
   do ispin=1,nspin
     write(file_name,'(a,i1,a)') 'rho_',ispin,'.cube'
     open(newunit=ocuberho(ispin),file=file_name)
     write(ocuberho(ispin),'(a)') 'cube file generated from MOLGW'
     write(ocuberho(ispin),'(a,i4)') 'density for spin ',ispin
     write(ocuberho(ispin),'(i6,3(f12.6,2x))') natom,xmin,ymin,zmin
     write(ocuberho(ispin),'(i6,3(f12.6,2x))') nx,dx,0.,0.
     write(ocuberho(ispin),'(i6,3(f12.6,2x))') ny,0.,dy,0.
     write(ocuberho(ispin),'(i6,3(f12.6,2x))') nz,0.,0.,dz
     do iatom=1,natom
       write(ocuberho(ispin),'(i6,4(2x,f12.6))') NINT(zatom(iatom)),0.0,x(:,iatom)
     enddo
   enddo
 endif

 do ix=1,nx
   rr(1) = xmin + (ix-1)*dx
   do iy=1,ny
     rr(2) = ymin + (iy-1)*dy
     do iz=1,nz
       rr(3) = zmin + (iz-1)*dz


       call calculate_basis_functions_r(basis,rr,basis_function_r)

       do ispin=1,nspin
         phi(istate1:istate2,ispin) = MATMUL( basis_function_r(:) , c_matrix(:,istate1:istate2,ispin) )
       enddo

       !
       ! check whether istate1:istate2 spans all the occupied states
       if( ALL( occupation(istate2+1,:) < completely_empty ) ) then
         do ispin=1,nspin
           write(ocuberho(ispin),'(50(e16.8,2x))') SUM( phi(:,ispin)**2 * occupation(istate1:istate2,ispin) )
         enddo
       endif


       do istate=istate1,istate2
         do ispin=1,nspin
           write(ocubefile(istate,ispin),'(50(e16.8,2x))') phi(istate,ispin)
         enddo
       enddo

     enddo
   enddo
 enddo

 deallocate(phi)

 do ispin=1,nspin
   do istate=istate1,istate2
     close(ocubefile(istate,ispin))
   enddo
   close(ocuberho(ispin))
 enddo

end subroutine plot_cube_wfn