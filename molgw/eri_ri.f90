!=========================================================================
function eri_ri(ibf,jbf,kbf,lbf)
  USE m_definitions
  USE m_eri
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_ri
!=====
 integer            :: index_ij,index_kl
 real(dp)           :: eri_1(1,1)
!=====
  LOGICAL :: negligible_basispair
  INTEGER :: index_pair

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_ri = 0.0_dp
 else
   index_ij = index_pair(ibf,jbf)
   index_kl = index_pair(kbf,lbf)
  
   eri_ri = DOT_PRODUCT( eri_3center(:,index_ij) , eri_3center(:,index_kl) )

   call xsum_auxil(eri_ri)

 endif

end function eri_ri
