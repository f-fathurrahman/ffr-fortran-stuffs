MODULE PSEUDOPOTENTIAL
implicit none
! number of species
integer                         :: N_species
! number of atoms per species
integer, allocatable            :: N_atom(:),z_atom(:)
integer                         :: N_atom_max
! positions
double precision,allocatable    :: Position(:,:,:)
integer                         :: N_orbitals
integer                         :: N_electrons,N_conduction
! max L in PS
integer,parameter               :: L_pp_max=4
! max number of radial points
integer,parameter               :: N_pp_max=1000
! pseudocore charge
double precision,allocatable    :: charge_pp(:)
! gaussian
double precision,parameter      :: beta=1.d0
! PP
double precision,allocatable    :: vion(:,:,:),p(:,:,:),r(:,:),clog(:)
integer,allocatable             :: l_loc(:),l_max(:),N_pp(:)
double precision,allocatable    :: rhops(:,:),vps(:,:)
complex*16      ,allocatable    :: fnl(:,:,:,:,:),ei1(:,:,:),ei2(:,:,:),ei3(:,:,:),eigr(:,:,:,:),sfac(:,:)
double precision,allocatable    :: wnl(:,:)
double precision,allocatable    :: pkg_a(:,:,:,:)



CONTAINS

subroutine input_atoms
  implicit none
  integer        :: i,j,k

  write(6,*)'?1'
  N_atom_max=0
  open(1,file='atom.inp')
  read(1,*)N_electrons,N_conduction
  read(1,*)n_species
  allocate(N_atom(N_species),charge_pp(N_species),z_atom(N_species))
  do i=1,n_species
    read(1,*)n_atom(i),z_atom(i)
    if(n_atom(i).gt.N_atom_max) N_atom_max=n_atom(i)
  end do
  allocate(Position(3,N_atom_max,N_species))
  do i=1,N_species
    do j=1,n_atom(i)
      read(1,*)(Position(k,j,i),k=1,3)
    end do
  end do
  N_orbitals=N_electrons/2+N_conduction
  if((-1)**(N_electrons/2).lt.0) N_orbitals=N_electrons/2+1+N_conduction
end subroutine input_atoms





subroutine read_pp
  implicit none
  integer          :: total_charge,i,l,j,is

  allocate (l_loc(N_species),l_max(N_species))
  allocate (vion(N_pp_max,N_species,L_pp_max))
  allocate (N_pp(N_species),clog(N_species))
  allocate (p(N_pp_max,N_species,3),r(N_pp_max,N_species))
  
  total_charge = 0
  i=0
  do is=1,n_species

    if(z_atom(is).eq.1) open(1,file='hpp.dat')
    if(z_atom(is).eq.3) open(1,file='lipp.dat')
    if(z_atom(is).eq.4) open(1,file='bepp.dat')
    if(z_atom(is).eq.5) open(1,file='bpp.dat')
    if(z_atom(is).eq.6) open(1,file='cpp.dat')
    if(z_atom(is).eq.7) open(1,file='npp.dat')
    if(z_atom(is).eq.8) open(1,file='opp.dat')
    if(z_atom(is).eq.11) open(1,file='napp.dat')
    if(z_atom(is).eq.13) open(1,file='alpp.dat')
    if(z_atom(is).eq.14) open(1,file='sipp.dat')
    if(z_atom(is).eq.16) open(1,file='spp.dat')
    if(z_atom(is).eq.28) open(1,file='nipp.dat')
    if(z_atom(is).eq.31) open(1,file='gapp.dat')
    if(z_atom(is).eq.33) open(1,file='aspp.dat')

    read(1,*) charge_pp(is), l_max(is), l_loc(is)

    total_charge = total_charge + n_atom(is) * charge_pp(is)

    do l=1,l_max(is)
      read(1,*) N_pp(is),clog(is),(i,r(j,is),p(j,is,l),vion(j,is,l),j=1,N_pp(is))
    enddo
    clog(is) = log(clog(is))
    close(1)
  enddo

end subroutine read_pp

subroutine initialize_pp
  USE GVECTOR
  USE PW_SMALL
  implicit none

  allocate (rhops(N_species,N_G_vector_max))
  allocate (vps(N_species,N_G_vector_max))
  allocate (fnl(N_G_wf_max,N_K_points,N_species,N_atom_max,(L_PP_max+1)**2))
  allocate (wnl(N_species,(L_PP_max+1)**2))
  allocate (pkg_a((L_PP_max+1)**2,N_G_K_vector_max,N_species,N_K_points))
  allocate (ei1(-(N_L(1)+1)/2:(N_L(1)+1)/2,N_atom_max,N_species))
  allocate (ei2(-(N_L(2)+1)/2:(N_L(2)+1)/2,N_atom_max,N_species))
  allocate (ei3(-(N_L(3)+1)/2:(N_L(3)+1)/2,N_atom_max,N_species))
  allocate (eigr(N_G_K_vector_max,N_atom_max,N_species,N_k_points))
  allocate (sfac(N_species,N_G_vector_max))

  call read_pp

end subroutine initialize_pp

END   MODULE PSEUDOPOTENTIAL
