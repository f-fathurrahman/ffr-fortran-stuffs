!=========================================================================
subroutine prepare_basis_functions_gradr(basis,batch_size)
  USE m_definitions
  USE m_basis_set
  USE m_dft_grid
 implicit NONE
 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: batch_size
!=====
 integer                    :: igrid
!=====

 write(stdout,*) 'Precalculate the gradients on N grid points',ngrid_stored
 if( batch_size /= 1 ) then
   write(stdout,'(3x,a,i5,a,i4)') 'which corresponds to ',ngrid_stored/batch_size,' batches of size ',batch_size
 endif
 call clean_allocate('basis grad ftns on grid',bfgr,basis%nbf,ngrid_stored,3)


 do igrid=1,ngrid_stored,batch_size
   call calculate_basis_functions_gradr_batch(basis,batch_size,rr_grid(:,igrid:igrid+batch_size-1),bfgr(:,igrid:igrid+batch_size-1,:))
 enddo

end subroutine prepare_basis_functions_gradr
