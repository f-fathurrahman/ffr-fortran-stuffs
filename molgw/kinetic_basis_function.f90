!=========================================================================
subroutine kinetic_basis_function(bf1,bf2,kinetic)
  use m_definitions, only : dp
  use m_basis_set
 implicit none
 type(basis_function),intent(in) :: bf1,bf2
 real(dp),intent(out)            :: kinetic
!=====
 integer                         :: ig,jg
 real(dp)                        :: kinetic_one_gaussian
!=====

 kinetic=0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call kinetic_recurrence(bf1%g(ig),bf2%g(jg),kinetic_one_gaussian)
     kinetic = kinetic + kinetic_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine kinetic_basis_function