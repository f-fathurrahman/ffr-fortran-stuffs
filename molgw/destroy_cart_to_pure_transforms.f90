!=========================================================================
subroutine destroy_cart_to_pure_transforms()
  use m_definitions
  use m_cart_to_pure
 implicit none

!=====
 integer :: il
!=====

 do il=0,MOLGW_LMAX
   deallocate(cart_to_pure_norm(il,CARTG)%matrix)
   deallocate(cart_to_pure_norm(il,PUREG)%matrix)
 enddo

end subroutine destroy_cart_to_pure_transforms
