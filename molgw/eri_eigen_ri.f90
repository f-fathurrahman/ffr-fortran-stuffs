!=========================================================================
function eri_eigen_ri(istate,jstate,ijspin,kstate,lstate,klspin)
  use m_definitions
  use m_eri
  use m_eri_ao_mo
 implicit none
 integer,intent(in) :: ijspin,klspin
 integer,intent(in) :: istate,jstate,kstate,lstate
 real(dp)           :: eri_eigen_ri
!=====

 eri_eigen_ri = DOT_PRODUCT( eri_3center_eigen(:,istate,jstate,ijspin) , eri_3center_eigen(:,kstate,lstate,klspin) )

 call xsum_auxil(eri_eigen_ri)

end function eri_eigen_ri
