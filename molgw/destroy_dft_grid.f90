!=========================================================================
subroutine destroy_dft_grid()
  use m_dft_grid
 implicit none

 deallocate(rr_grid)
 deallocate(w_grid)
 if( ALLOCATED(bf_rad2) ) then
   deallocate(bf_rad2)
 endif

 if( ALLOCATED(bfr) ) then
   call clean_deallocate('basis ftns on grid',bfr)
 endif
 if( ALLOCATED(bfgr) ) then
   call clean_deallocate('basis grad ftns on grid',bfgr)
 endif
 call destroy_dft_grid_distribution()
 
end subroutine destroy_dft_grid