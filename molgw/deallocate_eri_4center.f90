!=========================================================================
subroutine deallocate_eri_4center()
  USE m_eri, ONLY : eri_4center
 implicit none
!=====

 if(ALLOCATED(eri_4center)) then
   call clean_deallocate('4-center integrals',eri_4center)
 endif

end subroutine deallocate_eri_4center
