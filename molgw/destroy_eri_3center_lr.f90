!=================================================================
subroutine destroy_eri_3center_lr()
  USE m_eri
 implicit none
!=====

 if(ALLOCATED(iproc_ibf_auxil_lr)) then
   deallocate(iproc_ibf_auxil_lr)
 endif
 if(ALLOCATED(ibf_auxil_g_lr)) then
   deallocate(ibf_auxil_g_lr)
 endif
 if(ALLOCATED(ibf_auxil_l_lr)) then
   deallocate(ibf_auxil_l_lr)
 endif
 if(ALLOCATED(eri_3center_lr)) then
   call clean_deallocate('3-center LR integrals',eri_3center_lr)
 endif

end subroutine destroy_eri_3center_lr
