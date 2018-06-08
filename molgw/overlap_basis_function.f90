!=========================================================================
SUBROUTINE overlap_basis_function(bf1,bf2,overlap)
  USE m_definitions, only : dp
  USE m_basis_set
  IMPLICIT NONE 
  type(basis_function),intent(in) :: bf1,bf2
  real(dp),intent(out)            :: overlap
  !=====
  INTEGER :: ig,jg
  REAL(dp) :: overlap_one_gaussian
  !=====

 overlap=0.0_dp
 DO ig=1,bf1%ngaussian
   DO jg=1,bf2%ngaussian
     call overlap_recurrence(bf1%g(ig),bf2%g(jg),overlap_one_gaussian)
     overlap = overlap + overlap_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   ENDDO 
 ENDDO 


END SUBROUTINE overlap_basis_function


