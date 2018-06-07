!=========================================================================
function index_eri(ibf,jbf,kbf,lbf)
  USE m_eri
 implicit none

 integer,intent(in) :: ibf,jbf,kbf,lbf
 integer            :: index_eri
!=====
 integer            :: klmin,ijmax
 integer            :: index_ij,index_kl
!===== 
  INTEGER :: index_pair

 index_ij = index_pair(ibf,jbf)
 index_kl = index_pair(kbf,lbf)

 ijmax=MAX(index_ij,index_kl)
 klmin=MIN(index_ij,index_kl)

 index_eri = (klmin-1)*npair - (klmin-1)*(klmin-2)/2 + ijmax-klmin+1


end function index_eri
