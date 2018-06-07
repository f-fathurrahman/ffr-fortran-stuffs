!=========================================================================
subroutine basis_function_dipole(bf1,bf2,dipole)
  use m_definitions, only : dp
  use m_basis_set
 implicit none
 type(basis_function),intent(in)  :: bf1,bf2
 real(dp),intent(out)             :: dipole(3)
!=====
 type(basis_function)             :: bftmp
 real(dp)                         :: dipole_tmp
 logical,parameter                :: normalized=.FALSE.
 integer                          :: fake_shell=1
 integer                          :: fake_index=1
!=====

 ! 
 ! Calculate < phi_1 | r | phi_2 >
 ! using r = ( r - B ) + B
 !

 ! first set up | (x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(1) = dipole_tmp
 ! first set up | Bx phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bx phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(1) = dipole(1) + dipole_tmp * bf2%x0(1)

 ! first set up | (y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(2) = dipole_tmp
 ! first set up | By phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | By phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(2) = dipole(2) + dipole_tmp * bf2%x0(2)

 ! first set up | (z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(3) = dipole_tmp
 ! first set up | Bz phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bz phi2 >
 call overlap_basis_function(bf1,bftmp,dipole_tmp)
 dipole(3) = dipole(3) + dipole_tmp * bf2%x0(3)


end subroutine basis_function_dipole