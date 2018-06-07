!=========================================================================
subroutine deallocate_eri_4center_lr()
  USE m_eri, ONLY : eri_4center_lr
 implicit none
!=====

 if(ALLOCATED(eri_4center_lr)) then
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif

end subroutine deallocate_eri_4center_lr

