MODULE PW
  implicit none

! density
  complex*16,allocatable          :: store_density(:,:,:)
!  double precision,allocatable    :: density(:,:,:)
  complex*16,allocatable          :: HPsi(:)

  complex*8,allocatable           :: wave_function_c(:,:,:)
  double precision,allocatable    :: Occupation(:,:)
  double precision,allocatable    :: eigenvalue(:,:)
!
  integer                         :: N_mesh


  double precision                :: E_total
  double precision                :: E_kinetic
  double precision                :: E_Hartree
  double precision                :: E_exchange
  double precision                :: E_non_local
  double precision                :: E_self
  double precision                :: E_eigen
  double precision                :: E_pseudo
  double precision                :: E_es
  double precision                :: E_Gauss
  double precision                :: E_Fermi
  double precision                :: dens_mix

  double precision,parameter      :: rydberg=13.6058d0,pi=3.1415926535897d0
  double precision,parameter      :: kt=0.1d0

  integer,parameter               :: N_init=1

  logical                         :: Band_Structure

  CONTAINS

  subroutine initialize
  USE GVECTOR
  USE PW_SMALL
  USE PSEUDOPOTENTIAL
  USE FFT_DATA
  implicit none
  integer                         :: i,ik

  call input
  call input_atoms

  call initialize_lattice
  call initialize_symmetry

  call generate_G_vectors
  call pw_small_init
  call initialize_pp  

  N_mesh=(N_L(1)+fftinc1)*N_L(2)*N_L(3)
  allocate(store_density((N_L(1)+fftinc1),N_L(2),N_L(3)))
!  allocate(density((N_L(1)+fftinc1),N_L(2),N_L(3)))
  allocate(HPsi(N_G_K_vector_max))
  allocate (wave_function_c(N_G_K_vector_max,N_orbitals,N_K_points))
  allocate (eigenvalue(N_orbitals,N_K_points))
  allocate (occupation(N_orbitals,N_K_points))

  write(16,*)N_orbitals
  
  occupation=0.d0
  do i=1,N_electrons/2
    do ik=1,N_K_points
      occupation(i,ik)=2.0d0
    enddo
  enddo
  if((-1)**(N_electrons/2).lt.0) then
    do ik=1,N_K_points
      occupation(N_electrons/2+1,ik)=1.0d0
    enddo
  endif


  wave_function_c=(0.0,0.0)
  E_kinetic=0.d0
  dens_mix=0.95d0

  end subroutine initialize


  subroutine input
  USE GVECTOR
  USE PW_SMALL
  implicit none
  open(1,file='pw.inp')
    read(1,*) E_cut,E_cut_small
    read(1,*) BAND_STRUCTURE
  close(1)
    
  
  end subroutine input



END  MODULE PW
