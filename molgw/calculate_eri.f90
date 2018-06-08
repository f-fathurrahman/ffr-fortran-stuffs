!=========================================================================
subroutine calculate_eri(print_eri_,basis,rcut)
  USE m_definitions
  USE m_eri
  USE m_eri_calculate
  use m_timing
  use m_basis_set
 implicit none
 logical,intent(in)           :: print_eri_
 type(basis_set),intent(in)   :: basis
 real(dp),intent(in)          :: rcut
!=====
  logical :: read_eri

 call start_clock(timing_eri_4center)

 write(stdout,'(/,a,i12)') ' Number of integrals to be stored: ',nsize

 if( rcut < 1.0e-12_dp ) then
   call clean_allocate('4-center integrals',eri_4center,nsize)
   eri_4center(:) = 0.0_dp
 else
   call clean_allocate('4-center LR integrals',eri_4center_lr,nsize)
   eri_4center_lr(:) = 0.0_dp
 endif

 if( .NOT. read_eri(rcut) ) then
   call calculate_eri_4center(basis,rcut)
   if( print_eri_ ) then
     call dump_out_eri(rcut)
   endif
 endif

 call stop_clock(timing_eri_4center)

end subroutine calculate_eri
