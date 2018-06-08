!=========================================================================
function eri_lr(ibf,jbf,kbf,lbf)
  USE m_definitions
 implicit none
 integer,intent(in) :: ibf,jbf,kbf,lbf
 real(dp)           :: eri_lr
!=====
  LOGICAL :: negligible_basispair
  INTEGER :: index_eri
  REAL(dp) :: eri_4center_lr

 if( negligible_basispair(ibf,jbf) .OR. negligible_basispair(kbf,lbf) ) then
   eri_lr = 0.0_dp
 else
   eri_lr = eri_4center_lr(index_eri(ibf,jbf,kbf,lbf))
 endif

end function eri_lr
