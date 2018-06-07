!=================================================================
subroutine destroy_eri_3center_eigen()
  use m_eri_ao_mo
  implicit none
!=====

 write(stdout,'(/,a)') ' Destroy 3-center integrals on eigenstates'
 call clean_deallocate('3-center MO integrals',eri_3center_eigen)

end subroutine destroy_eri_3center_eigen
