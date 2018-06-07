!=========================================================================
subroutine print_basis_function(bf)
  use m_basis_set
 implicit none
 type(basis_function),intent(in) :: bf
!=====
 integer :: ig
!=====

 write(stdout,*)
 write(stdout,*) '======  print out a basis function ======'
 write(stdout,'(a30,2x,1(1x,i3))')           'contraction of N gaussians',bf%ngaussian
 write(stdout,'(a30,5x,a1)')                'orbital momentum',bf%amc
 write(stdout,'(a30,1x,3(f12.6,2x))')        'centered in',bf%x0(:)
 do ig=1,bf%ngaussian
   write(stdout,'(a30,2x,1x,i3,2x,f12.6)')   'coefficient',ig,bf%coeff(ig)
 enddo
 write(stdout,*)
 do ig=1,bf%ngaussian
   call print_gaussian(bf%g(ig))
 enddo
 write(stdout,*) '====== end of basis function ======'
 write(stdout,*)

end subroutine print_basis_function