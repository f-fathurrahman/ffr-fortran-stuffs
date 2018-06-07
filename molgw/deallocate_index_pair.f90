!=========================================================================
subroutine deallocate_index_pair()
  USE m_eri, ONLY : index_pair_1d
 implicit none

!=====

 if(ALLOCATED(index_pair_1d)) call clean_deallocate('index pair',index_pair_1d)

end subroutine deallocate_index_pair
