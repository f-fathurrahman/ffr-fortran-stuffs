!=========================================================================
subroutine dft_approximate_vhxc(print_matrix_,basis,vhxc_ij)
  use m_definitions
  use m_basis_set
  use m_timing
  use m_mpi
 use m_dft_grid
 use m_eri_calculate
 use m_atoms
  use m_inputparam, only : nspin
 implicit none

 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(out)       :: vhxc_ij(basis%nbf,basis%nbf)
!=====
 integer,parameter    :: BATCH_SIZE=64
 real(dp),allocatable :: weight_batch(:)
 real(dp),allocatable :: basis_function_r_batch(:,:)
 real(dp),allocatable :: exc_batch(:)
 real(dp),allocatable :: rhor_batch(:)
 real(dp),allocatable :: vrho_batch(:)
 integer              :: idft_xc
 integer              :: igrid,ibf,jbf,ispin
 integer              :: igrid_start,igrid_end
 integer              :: ir,nr
 real(dp)             :: normalization
 real(dp)             :: basis_function_r(basis%nbf)
 real(dp)             :: rhor
 real(dp)             :: vxc,excr,exc
 real(dp)             :: vsigma(2*nspin-1)
 real(dp)             :: vhgau(basis%nbf,basis%nbf)
 integer              :: iatom,igau,ngau
 real(dp),allocatable :: alpha(:),coeff(:)
!=====

 call start_clock(timing_approx_ham)

 vhxc_ij(:,:) = 0.0_dp

 write(stdout,'(/,a)') ' Calculate approximate HXC potential with a superposition of atomic densities'

 do iatom=1,natom
   if( rank_grid /= MODULO(iatom,nproc_grid) ) cycle

   ngau = 4
   allocate(alpha(ngau),coeff(ngau))
   call element_atomicdensity(zatom(iatom),basis_element(iatom),coeff,alpha)


   call calculate_eri_approximate_hartree(basis,basis%nbf,basis%nbf,x(:,iatom),ngau,coeff,alpha,vhgau)
   vhxc_ij(:,:) = vhxc_ij(:,:) + vhgau(:,:) 

   deallocate(alpha,coeff)
 enddo

 write(stdout,*) 'Simple LDA functional on a coarse grid'
 !
 ! Create a temporary grid with low quality
 ! This grid is to be destroyed at the end of the present subroutine
 call init_dft_grid(basis,low,.FALSE.,.FALSE.,BATCH_SIZE)


 normalization = 0.0_dp
 exc           = 0.0_dp
 do igrid_start=1,ngrid,BATCH_SIZE
   igrid_end = MIN(ngrid,igrid_start+batch_size-1)
   nr = igrid_end - igrid_start + 1

   allocate(weight_batch(nr))
   allocate(basis_function_r_batch(basis%nbf,nr))
   allocate(exc_batch(nr))
   allocate(rhor_batch(nr))
   allocate(vrho_batch(nr))

   weight_batch(:) = w_grid(igrid_start:igrid_end)

   !
   ! Get all the functions and gradients at point rr
   call get_basis_functions_r_batch(basis,igrid_start,nr,basis_function_r_batch)

   !
   ! Calculate the density at points r
   do ir=1,nr
     igrid = igrid_start + ir - 1
     call setup_atomic_density(rr_grid(:,igrid),rhor_batch(ir))
   enddo

   !
   ! Normalization
   normalization = normalization + SUM( rhor_batch(:) * weight_batch(:) )

   call teter_lda_vxc_exc(nr,rhor_batch,vrho_batch,exc_batch)

   !
   ! XC energy
   exc = exc + SUM( exc_batch(:) * weight_batch(:) * rhor_batch(:) )

   !
   ! HXC
   ! Hopefully, vrho is always negative
   do ir=1,nr
     basis_function_r_batch(:,ir) = SQRT( -weight_batch(ir) * vrho_batch(ir) ) * basis_function_r_batch(:,ir)
   enddo
   call DSYRK('L','N',basis%nbf,nr,-1.0d0,basis_function_r_batch,basis%nbf,1.0d0,vhxc_ij,basis%nbf)

!   call DSYR('L',basis%nbf,weight*vxc,basis_function_r,1,vhxc_ij,basis%nbf)

   deallocate(weight_batch)
   deallocate(basis_function_r_batch)
   deallocate(exc_batch)
   deallocate(rhor_batch)
   deallocate(vrho_batch)

 enddo ! loop on the grid point
 !
 ! Sum up the contributions from all procs only if needed
 call xsum_grid(normalization)
 call xsum_grid(exc)
 call xsum_grid(vhxc_ij)

 ! Symmetrize vhxc
 do ibf=1,basis%nbf
   do jbf=ibf+1,basis%nbf
     vhxc_ij(ibf,jbf) = vhxc_ij(jbf,ibf)
   enddo
 enddo

 write(stdout,'(/,a,2(2x,f12.6))') ' Number of electrons:',normalization
 write(stdout,  '(a,2(2x,f12.6))') '      XC energy (Ha):',exc

 !
 ! Temporary grid destroyed
 call destroy_dft_grid()

 call stop_clock(timing_approx_ham)

end subroutine dft_approximate_vhxc
