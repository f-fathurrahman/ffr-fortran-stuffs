!=========================================================================
FUNCTION negligible_basispair(ibf,jbf)
  USE m_eri
 IMPLICIT NONE
 INTEGER, INTENT(in) :: ibf,jbf
 LOGICAL :: negligible_basispair
!=====
 INTEGER  :: ishell,jshell
!=====


 ishell = shell_bf(ibf)
 jshell = shell_bf(jbf)

 negligible_basispair = negligible_shellpair(ishell,jshell)


end function negligible_basispair
