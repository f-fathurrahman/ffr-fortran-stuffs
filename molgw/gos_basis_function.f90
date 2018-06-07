!=========================================================================
subroutine gos_basis_function(bf1,bf2,qvec,gos_bf1bf2)
  use m_definitions, only : dp
  use m_basis_set
 implicit none
 type(basis_function),intent(in)  :: bf1,bf2
 real(dp),intent(in)              :: qvec(3)
 complex(dp),intent(out)          :: gos_bf1bf2
!=====
 integer                          :: ig,jg
 complex(dp)                      :: gos_one_gaussian
!=====

 gos_bf1bf2 = 0.0_dp
 do ig=1,bf1%ngaussian
   do jg=1,bf2%ngaussian
     call evaluate_gos(bf1%g(ig),bf2%g(jg),qvec,gos_one_gaussian)
     gos_bf1bf2 = gos_bf1bf2 + gos_one_gaussian * bf1%coeff(ig) * bf2%coeff(jg)
   enddo
 enddo


end subroutine gos_basis_function
