!=========================================================================
subroutine deallocate_eri()
  USE m_eri
 implicit none

 integer :: ishell
!=====

 if(ALLOCATED(eri_4center)) then
   call clean_deallocate('4-center integrals',eri_4center)
 endif
 if(ALLOCATED(eri_4center_lr)) then
   call clean_deallocate('4-center LR integrals',eri_4center_lr)
 endif
 if(ALLOCATED(negligible_shellpair))   deallocate(negligible_shellpair)
 if(ALLOCATED(index_shellpair))        deallocate(index_shellpair)
 if(ALLOCATED(index_pair_1d)) call clean_deallocate('index pair',index_pair_1d)
 if(ALLOCATED(index_basis))   call clean_deallocate('index basis',index_basis)

 if(ALLOCATED(shell_bf))              deallocate(shell_bf)


end subroutine deallocate_eri
