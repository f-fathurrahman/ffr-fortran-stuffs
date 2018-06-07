!=========================================================================
subroutine overlap_basis_function(bf1,bf2,overlap)
  use m_definitions, only : dp
  use m_basis_set
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(out)            :: overlap
!=====
 integer                         :: ig,jg
 real(dp)                        :: overlap_one_gaussian
!=====

 overlap=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call overlap_recurrence(bf1%g(ig),bf2%g(jg),overlap_one_gaussian)
     overlap = overlap + overlap_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine overlap_basis_function