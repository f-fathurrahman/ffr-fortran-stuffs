!=========================================================================
subroutine prepare_basis_functions_r(basis,batch_size)
 implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: batch_size
!=====
 integer                    :: igrid
!=====

 write(stdout,'(1x,a,i7)') 'Precalculate the functions on N grid points ',ngrid_stored
 if( batch_size /= 1 ) then
   write(stdout,'(3x,a,i5,a,i4)') 'which corresponds to ',ngrid_stored/batch_size,' batches of size ',batch_size
 endif
 call clean_allocate('basis ftns on grid',bfr,basis%nbf,ngrid_stored)

 do igrid=1,ngrid_stored,batch_size
   call calculate_basis_functions_r_batch(basis,batch_size,rr_grid(:,igrid:igrid+batch_size-1),bfr(:,igrid:igrid+batch_size-1))
 enddo

end subroutine prepare_basis_functions_r
