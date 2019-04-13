MODULE GVECTOR
implicit none
! Lattice dimension
  integer                      :: N_L(3)
! cut off in G_space
  double precision             :: E_cut,G_cut
! max number of G_vectors
  integer                      :: N_G_vector_max
! G_vectors (must be multiplied by the B_Lattice vector
  integer,allocatable          :: G_vector(:,:)
  double precision,allocatable :: G_vector_length(:)
! K points
  integer                      :: N_K_points
  double precision,allocatable :: K_point(:,:),W_K_point(:)
! number of G vectors for a given K
  integer,allocatable          :: N_G_vector(:)
! max number of G_vectors for a given K
  integer                      :: N_G_K_vector_max
! index pointing to the position of the G vector 
  integer,allocatable          :: G_index(:,:)
! (G+k)**2
  double precision,allocatable :: Gplusk(:,:)
! fft index
  integer,allocatable          :: fft_ind(:,:),ind3f1(:,:),ind3f2(:,:),ind2f1(:),ind2f2(:)    
! Lattice  vector
  double precision             :: Lattice_vector(3,3)
  double precision             :: Lattice_constant
! Lattice vector length            
  double precision             :: Lattice_vector_length(3)
!
! Reciprocal space lattice
!
! Lattice  vector
  double precision             :: R_Lattice_vector(3,3),Volume
! Lattice vector length            
  double precision             :: R_Lattice_vector_length(3)
! symmetry
  integer                      :: N_sym,sym_mat(3,3,48)
!
  double precision             :: tpiba,tpiba2
  
  CONTAINS

  subroutine initialize_lattice
    implicit none
    integer         :: i,j
    real*8          :: g(3),a(3,3),b(3,3),pi
    pi=4.d0*atan(1.d0)


  open(1,file='lattice.inp')
    read(1,*)Lattice_constant
    do i=1,3
      read(1,*)(Lattice_vector(j,i),j=1,3)
    end do
  close(1)


  open(1,file='kpoints.inp')
    read(1,*)N_k_points
    allocate(K_point(3,N_K_Points),W_K_point(N_K_points))
    do i=1,N_k_points
      read(1,*)(K_point(j,i),j=1,3),W_K_point(i)
    end do
  close(1)
  
!
! rescale lattice vectors
!  Lattice_vector=Lattice_vector*Lattice_constant
  Lattice_vector=Lattice_vector
! length of lattice vectors
  do i=1,3
    Lattice_vector_length(i) = SQRT(dot_product(Lattice_vector(:,i),Lattice_vector(:,i)))
  end do
!
! reciprocal lattice
!
  a=Lattice_vector
  call inv_r(a,3,b)
  R_lattice_vector=Transpose(b)*Lattice_constant
! length of lattice vectors
  do i=1,3
    R_Lattice_vector_length(i) = SQRT(dot_product(R_Lattice_vector(:,i),R_Lattice_vector(:,i)))
  end do
  g(1)=Lattice_vector(2,2)*Lattice_vector(3,3)-Lattice_vector(3,2)*Lattice_vector(2,3)
  g(2)=Lattice_vector(3,2)*Lattice_vector(1,3)-Lattice_vector(1,2)*Lattice_vector(3,3)
  g(3)=Lattice_vector(1,2)*Lattice_vector(2,3)-Lattice_vector(2,2)*Lattice_vector(1,3)
  Volume=sum(Lattice_vector(:,1)*g(:))
  N_L=2*int(sqrt(e_cut)*Lattice_vector_length/Pi+1)
  write(33,*)R_Lattice_vector
  write(33,*)'N_L',N_L
  G_cut=e_cut*Lattice_constant**2/(2.d0*Pi)**2
  write(33,*)G_cut


  tpiba=2.0*pi/Lattice_constant
  tpiba2=tpiba**2


  end subroutine initialize_lattice

subroutine initialize_symmetry
  implicit none
  integer        :: i,j,k,i1,i2,i3
  integer        :: s(3,3,48)
  open(1,file='symmetry.inp')
    read(1,*)N_Sym
    do i=1,N_sym
      read(1,*)
      do k=1,3
        read(1,*)(s(j,k,i),j=1,3)
      end do
    end do
  close(1)
  
   do i=1,N_sym
     do i1=1,3  
       do i2=1,3 
         if(mod(s(i1,i2,i)*N_L(i2),N_L(i1)).ne.0) write(6,*)'wrong symmetry'
         sym_mat(i1,i2,i)=real(s(i1,i2,i)*N_L(i2))/real(N_L(i1))
       enddo
     enddo
   enddo
   do i=1,N_k_points
     W_K_point(i)=W_K_point(i)/dfloat(N_sym)
   end do

end subroutine initialize_symmetry

subroutine Generate_G_vectors
  implicit none
  integer                      :: i1,i2,i3,n1,n2,n3,ng,ii,i,j,ig(3),k
  integer,allocatable          :: igv(:,:),inc(:)
  double precision,allocatable :: gv(:)
  double precision             :: g(3),gg
  logical,allocatable          :: ic(:,:)
  logical                      :: con
  integer,external             :: iflip_inv,iflip
    
!   two cycles: first count and allocate then store 
  do ii=1,2
!   upper half: count the g vectors with g(1).ge.0
    ng=1
    do i3=0,N_L(3)/2
      n2=-N_L(2)/2
      if(i3.eq.0) n2=0
      do i2=n2,N_L(2)/2
        n1=-N_L(1)/2
        if (i3.eq.0 .and.i2.eq.0) n1=1
        do i1=n1,N_L(1)/2
          g(:)=i1*R_Lattice_vector(:,1)+i2*R_Lattice_vector(:,2)+i3*R_Lattice_vector(:,3)
          gg=sum(g(:)**2)
          if(gg.le.4.d0*G_cut) then
            ng=ng+1          
            if(ii.eq.2) then
              gv(ng)=gg
              igv(1,ng)=i1
              igv(2,ng)=i2
              igv(3,ng)=i3
            endif
          endif
        end do
      end do
    end do
!   upperhalf +lowerhalf+ g=0    
    N_G_vector_max=2*ng-1
!   now allocate arrays to store this
    if(ii.eq.1) then
      allocate(gv(ng),igv(3,ng),inc(ng))
!   G=0 component
      gv(1)=0.d0
      igv(:,1)=0
    endif    
  end do
!
! reorder with increasing magnitude
!
  call ordering(ng,gv,inc)
!  call indexx(ng,gv,inc)
  allocate(G_vector(3,N_G_vector_max),G_vector_length(N_G_vector_max))
! G=0
  G_vector(:,1)=1
  G_vector_length(1)=0.d0
  ii=2
  do i=2,ng
    do k=1,3
      G_vector(k,ii+1)=iflip_inv(-igv(k,inc(i)),N_L(k))
      G_vector(k,ii)=iflip_inv(igv(k,inc(i)),N_L(k))
    end do
    G_vector_length(ii)=gv(inc(i))
    G_vector_length(ii+1)=gv(inc(i))
    ii=ii+2
  end do
!
! G vectors for a given K point
!
  allocate(ic(N_L(1),N_L(2)))
  allocate(N_G_vector(N_K_points),G_index(N_G_vector_max,N_K_points))
  allocate(GplusK(N_G_vector_max,N_K_points),fft_ind(N_G_vector_max,N_K_points))
  allocate(ind3f1(N_L(1),N_K_points),ind3f2(N_L(1),N_K_points))
  allocate(ind2f1(N_K_points),ind2f2(N_K_points))

  N_G_K_vector_max=0
  
  do k=1,N_K_Points
    ic=.true.
    ii=0
    do i=1,N_G_vector_max
      ig(1)=iflip(G_vector(1,i),N_L(1))
      ig(2)=iflip(G_vector(2,i),N_L(2))
      ig(3)=iflip(G_vector(3,i),N_L(3))
      g(:) = ig(1)*R_Lattice_vector(:,1)+ig(2)*R_Lattice_vector(:,2)+ig(3)*R_Lattice_vector(:,3)+K_point(:,k)
      gg=sum(g(:)**2)
      if(gg.le.G_cut) then
        ii=ii+1
        G_index(ii,k)=i 
        Gplusk(ii,k)=gg
        fft_ind(ii,k)=G_vector(1,i)+(N_L(1)+1)*((G_vector(2,i)-1)+N_L(2)*(G_vector(3,i)-1))
        ic(G_vector(1,i),G_vector(2,i))=.false.
      endif
    end do
    N_G_vector(k)=ii
    if(N_G_vector(k).gt.N_G_K_vector_max) N_G_K_vector_max=N_G_vector(k)
    
!
!   index limits for the FFT for the wf and potential
!
    ind2f1(k)=1
    ind2f2(k)=N_L(1)+1
    do i=1,N_L(1)
      con=.true.
      ind3f1(i,k)=0
      ind3f2(i,k)=N_L(2)+1
      do j=N_L(2)/2,1,-1
        if ( .not.ic(i,j) .and. con) then
           ind3f1(i,k)=j
           con=.false.
         endif
      enddo
      con=.true.
      do j=N_L(2)/2,N_L(2)
        if (.not.ic(i,j).and.con) then
          ind3f2(i,k)=j
          con=.false.
        endif
      enddo
      if ((ind3f1(i,k).ne.0).and.(ind3f2(i,k).ne.N_L(2)+1).and.(ind2f2(k).eq.N_L(1)+1)) ind2f1(k)=i+1
      if ((ind3f1(i,k).eq.0).and.(ind3f2(i,k).eq.N_L(2)+1)) ind2f2(k)=i+1
    enddo
    
  end do
  
end subroutine Generate_G_vectors


END MODULE GVECTOR
