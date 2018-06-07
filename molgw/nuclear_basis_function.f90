!=========================================================================
subroutine nucleus_basis_function(bf1,bf2,zatom,x,nucleus_pot)
  use m_definitions, only : dp
  use m_basis_set
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(in)             :: zatom,x(3)
 real(dp),intent(out)            :: nucleus_pot
!=====
 integer                         :: ig,jg
 real(dp)                        :: nucleus_pot_one_gaussian
!=====

 nucleus_pot=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call nucleus_recurrence(zatom,x,bf1%g(ig),bf2%g(jg),nucleus_pot_one_gaussian)
     nucleus_pot = nucleus_pot + nucleus_pot_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine nucleus_basis_function