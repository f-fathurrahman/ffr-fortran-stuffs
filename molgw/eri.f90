!=========================================================================
function eri(ibf,jbf,kbf,lbf)
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri
!=====

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri = 0.0_dp
 else
   eri = eri_4center(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri
