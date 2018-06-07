!=========================================================================
subroutine basis_function_quadrupole(bf1,bf2,quad)
  use m_definitions, only : dp
  use m_basis_set
  implicit none
  type(basis_function),intent(in)  :: bf1,bf2
  real(dp),intent(out)             :: quad(3,3)
  !=====
  type(basis_function)             :: bftmp
  real(dp)                         :: quad_tmp
  logical,parameter                :: normalized=.FALSE.
  integer                          :: fake_shell=1
  integer                          :: fake_index=1
  !=====

 ! 
 ! Calculate < phi_1 | x y | phi_2 >
 ! using x y = ( x - Bx ) * ( y - By) + Bx * ( y - By ) + By * ( x - Bx ) + Bx * By
 !



 !
 !  terms x * y and y * x
 ! 

 ! first set up | (x-Bx)*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad_tmp

 ! first set up | Bx*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(1)

 ! first set up | By*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(2)

 ! Overlap < phi1 | Bx*By phi2 >
 call overlap_basis_function(bf1,bf2,quad_tmp)
 quad(1,2) = quad(1,2) + quad_tmp * bf2%x0(1) * bf2%x0(2)

 quad(2,1) = quad(1,2)


 !
 !  terms x * z and z * x
 ! 

 ! first set up | (x-Bx)*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad_tmp

 ! first set up | Bx*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bx*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(1)

 ! first set up | Bz*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bz*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(3)

 ! Overlap < phi1 | Bx*Bz phi2 >
 call overlap_basis_function(bf1,bf2,quad_tmp)
 quad(1,3) = quad(1,3) + quad_tmp * bf2%x0(1) * bf2%x0(3)

 quad(3,1) = quad(1,3)


 !
 !  terms y * z and z * y
 ! 

 ! first set up | (y-By)*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By)*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad_tmp

 ! first set up | By*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | By*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(2)

 ! first set up | Bz*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | Bz*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(3)

 ! Overlap < phi1 | By*Bz phi2 >
 call overlap_basis_function(bf1,bf2,quad_tmp)
 quad(2,3) = quad(2,3) + quad_tmp * bf2%x0(2) * bf2%x0(3)

 quad(3,2) = quad(2,3)




 !
 !  term x * x
 ! 

 ! first set up | (x-Bx)*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+2,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (x-Bx)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad_tmp

 ! first set up | 2Bx*(x-Bx) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx+1,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2Bx*(x-Bx) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad(1,1) + quad_tmp * 2.0_dp * bf2%x0(1)

 ! first set up | Bx**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bx**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(1,1) = quad(1,1) + quad_tmp * bf2%x0(1)**2


 !
 !  term y * y
 ! 

 ! first set up | (y-By)^2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+2,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (y-By)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad_tmp

 ! first set up | 2B*(y-By) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny+1,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2B*(y-By) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad(2,2) + quad_tmp * 2.0_dp * bf2%x0(2)

 ! first set up | By**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | By**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(2,2) = quad(2,2) + quad_tmp * bf2%x0(2)**2


 !
 !  term z * z
 ! 

 ! first set up | (z-Bz)^2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+2,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | (z-Bz)^2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad_tmp

 ! first set up | 2B*(z-Bz) phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz+1,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor
 ! then overlap < phi1 | 2B*(z-Bz) phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad(3,3) + quad_tmp * 2.0_dp * bf2%x0(3)

 ! first set up | Bz**2 phi_2 >
 call init_basis_function(normalized,bf2%ngaussian,bf2%nx,bf2%ny,bf2%nz,0,bf2%x0,bf2%g(:)%alpha,bf2%coeff,fake_shell,fake_index,bftmp)
 ! override the usual normalization
 bftmp%g(:)%norm_factor = bf2%g(:)%norm_factor 
 ! then overlap < phi1 | Bz**2 phi2 >
 call overlap_basis_function(bf1,bftmp,quad_tmp)
 quad(3,3) = quad(3,3) + quad_tmp * bf2%x0(3)**2

end subroutine basis_function_quadrupole
