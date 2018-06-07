!=================================================================
subroutine destroy_eri_3center_eigen_mixed()
  use m_eri_ao_mo
 implicit none
!=====

 write(stdout,'(/,a)') ' Destroy 3-center mixed integrals on eigenstates'
 call clean_deallocate('3-center mixed MO integrals',eri_3center_eigen_mixed)

end subroutine destroy_eri_3center_eigen_mixed
