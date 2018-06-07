!=================================================================
subroutine destroy_eri_3center()
  USE m_eri
 implicit none
!=====

 if(ALLOCATED(iproc_ibf_auxil)) then
   deallocate(iproc_ibf_auxil)
 endif
 if(ALLOCATED(ibf_auxil_g)) then
   deallocate(ibf_auxil_g)
 endif
 if(ALLOCATED(ibf_auxil_l)) then
   deallocate(ibf_auxil_l)
 endif
 if(ALLOCATED(eri_3center)) then
   call clean_deallocate('3-center integrals',eri_3center)
 endif

#ifdef SCASCA
 if(ALLOCATED(eri_3center_sca)) then
   call clean_deallocate('3-center integrals SCALAPACK',eri_3center_sca)
 endif
#endif

end subroutine destroy_eri_3center
