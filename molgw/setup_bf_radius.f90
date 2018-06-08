!=========================================================================
subroutine setup_bf_radius(basis)
  use m_definitions
  use m_basis_set
  use m_dft_grid
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer  :: igrid,ibf
 real(dp) :: basis_function_r(basis%nbf)
 integer  :: icalc,icalc_tot
!=====

 allocate(bf_rad2(basis%nbf))
 bf_rad2(:) = 0.0_dp
 do igrid=1,ngrid
   call get_basis_functions_r(basis,igrid,basis_function_r)
   do ibf=1,basis%nbf
     if( ABS(basis_function_r(ibf)) > TOL_BF ) then
       bf_rad2(ibf) = MAX( bf_rad2(ibf) , SUM( (rr_grid(:,igrid) - basis%bff(ibf)%x0(:))**2 ) )
     endif
   enddo
 enddo

 call xmax_grid(bf_rad2)

end subroutine setup_bf_radius
