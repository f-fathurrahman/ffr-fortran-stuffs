!=================================================================
subroutine destroy_eri_3center_eigen_lr()
  use m_eri_ao_mo
 implicit none
!=====

 write(stdout,'(/,a)') ' Destroy LR 3-center integrals on eigenstates'
 call clean_deallocate('LR 3-center MO integrals',eri_3center_eigen_lr)

end subroutine destroy_eri_3center_eigen_lr