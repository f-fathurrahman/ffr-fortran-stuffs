!=========================================================================
function index_pair(ibf,jbf)
  USE m_eri
 implicit none

 integer,intent(in) :: ibf,jbf
 integer            :: index_pair
!=====
 integer            :: ijmin,ijmax
!=====

 if( ibf == jbf ) then
   index_pair = ibf
 else
   ijmax=MAX(ibf,jbf)
   ijmin=MIN(ibf,jbf)

   index_pair = (ijmin-1) * (nbf_eri-1) - (ijmin-1) * (ijmin-2)/2     + ijmax - ijmin + nbf_eri

   index_pair = index_pair_1d(index_pair)
 endif


end function index_pair
