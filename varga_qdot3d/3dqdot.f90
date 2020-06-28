!----------------------------
SUBROUTINE poisson_solve_3d()
!----------------------------
  USE m_qd3d
  IMPLICIT NONE

  integer :: k1,k2,k3

  rho = -4.d0*pi*rho
  
  ! Boundary condition in the x direction
  do k2=1,N_L(2)
    do k3=1,N_L(3)
      rho(Lattice_inv(1,k2,k3))=rho(Lattice_inv(1,k2,k3))-V_X0(k2,k3)/grid_step(1)**2/E2
      rho(Lattice_inv(N_L(1),k2,k3))=rho(Lattice_inv(N_L(1),k2,k3))-V_XN(k2,k3)/grid_step(1)**2/E2
    enddo
  enddo
  
  ! Boundary condition in the y direction
  do k1=1,N_L(1)
    do k3=1,N_L(3)
      rho(Lattice_inv(k1,1,k3))=rho(Lattice_inv(k1,1,k3))-V_Y0(k1,k3)/grid_step(2)**2/E2
      rho(Lattice_inv(k1,N_L(2),k3))=rho(Lattice_inv(k1,N_L(2),k3))-V_YN(k1,k3)/grid_step(2)**2/E2
    enddo
  enddo
  
  ! Boundary condition in the z direction
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      rho(Lattice_inv(k1,k2,1))=rho(Lattice_inv(k1,k2,1))-V_Z0(k1,k2)/grid_step(3)**2/E2
      rho(Lattice_inv(k1,k2,N_L(3)))=rho(Lattice_inv(k1,k2,N_L(3)))-V_ZN(k1,k2)/grid_step(3)**2/E2
    enddo
  enddo

  call CoGr(rho,N_L_points)
  v_pot=E2*rho
END SUBROUTINE


!------------------------
SUBROUTINE init_lattice()
!------------------------
  USE m_qd3d
  IMPLICIT NONE  
  integer :: k1,k2,k3,num,i,k

  ! Setup the lattice bookkeeping
  num=0
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      do k3=1,N_L(3)
        num=num+1 
        Lattice(1,num)=k1
        Lattice(2,num)=k2
        Lattice(3,num)=k3
        Lattice_inv(k1,k2,k3)=num
        grid_point(1,num)=-0.5d0*(N_L(1)-1)*grid_step(1)+(k1-1)*grid_step(1)
        grid_point(2,num)=-0.5d0*(N_L(2)-1)*grid_step(2)+(k2-1)*grid_step(2)
        grid_point(3,num)=-0.5d0*(N_L(3)-1)*grid_step(3)+(k3-1)*grid_step(3)
      enddo
    enddo
  enddo
END SUBROUTINE
  
!----------------------------------
SUBROUTINE init_confining_potential
!----------------------------------
  USE m_qd3d
  IMPLICIT NONE

  ! Setup the confining potential
  integer :: i,k1,k2,k3
  REAL(8)  :: r0(3),x2
  
  r0=(/0.d0,0.d0,0.d0/)
 
  do i=1,N_L_points
    x2=(grid_point(1,i)-r0(1))**2+(grid_point(2,i)-r0(2))**2+(grid_point(3,i)-r0(3))**2
    V_ext(i)=0.5d0*omega0**2*x2!+10.d0*exp(-x2*10.d0)
  enddo
   
  k2=int(dfloat(N_L(2))/2.d0)
  k3=int(dfloat(N_L(3))/2.d0)
  do k1=1,N_L(1)
	  i=Lattice_inv(k1,k2,k3)
    write(666,*)grid_point(1,i),V_ext(i)
  enddo  
END SUBROUTINE


!-----------------------------
SUBROUTINE laplace_operator4()
!-----------------------------
  USE m_qd3d
  IMPLICIT NONE

  ! 4th order finite difference representation of the Laplacian
  REAL(8), parameter :: cN0=-205.d0/72.d0,cN1=8.d0/5.d0
  REAL(8), parameter :: cN2=-1.d0/5.d0,cN3=8.d0/315.d0, cN4=-1.d0/560.d0
  integer          :: i,i1,i2,i3
  REAL(8)           :: K_x,K_y,K_z

  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    wf(i1,i2,i3)=Phi(i)
  enddo

  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    K_x=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1+1,i2,i3)+wf(i1-1,i2,i3))+& 
             cN2*(wf(i1+2,i2,i3)+wf(i1-2,i2,i3))+& 
             cN3*(wf(i1+3,i2,i3)+wf(i1-3,i2,i3))+& 
             cN4*(wf(i1+4,i2,i3)+wf(i1-4,i2,i3)))/Grid_Step(1)**2
    K_y=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1,i2+1,i3)+wf(i1,i2-1,i3))+& 
             cN2*(wf(i1,i2+2,i3)+wf(i1,i2-2,i3))+& 
             cN3*(wf(i1,i2+3,i3)+wf(i1,i2-3,i3))+& 
             cN4*(wf(i1,i2+4,i3)+wf(i1,i2-4,i3)))/Grid_Step(2)**2
    K_z=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1,i2,i3+1)+wf(i1,i2,i3-1))+& 
             cN2*(wf(i1,i2,i3+2)+wf(i1,i2,i3-2))+& 
             cN3*(wf(i1,i2,i3+3)+wf(i1,i2,i3-3))+& 
             cN4*(wf(i1,i2,i3+4)+wf(i1,i2,i3-4)))/Grid_Step(3)**2
    L_Phi(i)=K_x+K_y+K_z
  enddo
END SUBROUTINE

!----------------------------
SUBROUTINE laplace_operator()
!----------------------------
  USE m_qd3d
  IMPLICIT NONE

  ! 1st order finite difference representation of the Laplacian
  integer :: i1,i2,i3,kpoint,i
  REAL(8)  :: k_x,K_y,K_z
  wf=0.d0    
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    wf(i1,i2,i3)=phi(i)
  end do
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    K_x=(-2.d0*wf(i1,i2,i3)+wf(i1+1,i2,i3)+wf(i1-1,i2,i3))/grid_step(1)**2
    K_y=(-2.d0*wf(i1,i2,i3)+wf(i1,i2+1,i3)+wf(i1,i2-1,i3))/grid_step(2)**2
    K_z=(-2.d0*wf(i1,i2,i3)+wf(i1,i2,i3+1)+wf(i1,i2,i3-1))/grid_step(3)**2
    L_phi(i)=K_x+K_y+K_z
  enddo
END SUBROUTINE laplace_operator

!--------------------
SUBROUTINE CoGr(b,NL)
!--------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Solution of the Poisson equation using the conjugate gradient method
	integer                          :: NL,iteration,i,j
	REAL(8)                           :: c0,c1,alfa,beta,rr,bn,con,b(NL)
	integer,parameter                :: N_iter=1500
	REAL(8),parameter                 :: eps=1.d-10
	REAL(8),dimension(:),allocatable  :: g,d,h,x

	allocate(g(NL),d(NL),h(NL),x(NL))

	bn=sqrt(dot_product(b,b))
	x=0.d0
	phi=x
	if(N_d==4) then
	  call laplace_operator4()
	else
	  call laplace_operator()
	endif
	g = b + L_phi
	d = g
	c0 = dot_product(g,g)

	do iteration=1,N_iter
		con=abs(sqrt(dot_product(g,g))/bn)
		if(con>eps) then
		  phi=d
		  	if(N_d==4) then
					call laplace_operator4
				else
					call laplace_operator
				endif
		  h=L_phi
		  alfa=-c0/dot_product(d,h)
		  x=x+alfa*d
		  g=g+alfa*h
		  c1=dot_product(g,g)
		  beta=c1/c0; c0=c1
		  d=g+beta*d
		endif
	enddo
	b=x
	if(con>eps) then
		write(*,*)'Poisson is not converged!'
	endif
	deallocate(g,d,x,h)
END SUBROUTINE

!------------------------
SUBROUTINE bc_multipole()
!------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Boundary conditions determined by multipole expansion 
  integer :: k1,k2,k3,i,l,lm
  REAL(8) :: xx,yy,zz,x,y,z,r,rho_lm((L_max+1)**2)
  ! Function
  REAL(8) :: Ylm, mp_pot

  rho_lm=0.d0
  do L=0,L_max
    do lm=L**2+1,(L+1)**2
      do i=1,N_L_points
         x = grid_point(1,i)
         y = grid_point(2,i)
         z = grid_point(3,i)
         r = sqrt(x*x+y*y+z*z) + small
         xx=x/r
         yy=y/r
         zz=z/r
         rho_lm(lm) = rho_lm(lm) + r**L * Ylm(xx,yy,zz,lm) * rho(i) * grid_volume
       enddo
    enddo
  enddo
  
  ! Boundary condition in the x direction
  do k2=1,N_L(2)
    do k3=1,N_L(3)
      y=grid_point(2,Lattice_inv(1,k2,k3))
      z=grid_point(3,Lattice_inv(1,k2,k3))
      x=grid_point(1,Lattice_inv(1,k2,k3))-grid_step(1)
      V_X0(k2,k3)=mp_pot(x,y,z,rho_lm)
      x=grid_point(1,Lattice_inv(N_L(1),k2,k3))+grid_step(1)
      V_XN(k2,k3)=mp_pot(x,y,z,rho_lm)
    enddo
  enddo
  
  ! Boundary condition in the y direction
  do k1=1,N_L(1)
    do k3=1,N_L(3)
      x=grid_point(1,Lattice_inv(k1,1,k3))
      z=grid_point(3,Lattice_inv(k1,1,k3))
      y=grid_point(2,Lattice_inv(k1,1,k3))-grid_step(2)
      V_Y0(k1,k3)=mp_pot(x,y,z,rho_lm)
      y=grid_point(2,Lattice_inv(k1,N_L(2),k3))+grid_step(2)
      V_YN(k1,k3)=mp_pot(x,y,z,rho_lm)
    enddo
  enddo
  
  ! Boundary condition in the z direction
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      x=grid_point(1,Lattice_inv(k1,k2,1))
      y=grid_point(2,Lattice_inv(k1,k2,1))
      z=grid_point(3,Lattice_inv(k1,k2,1))-grid_step(3)
      V_Z0(k1,k2)=mp_pot(x,y,z,rho_lm)
      z=grid_point(3,Lattice_inv(k1,k2,N_L(3)))+grid_step(3)
      V_ZN(k1,k2)=mp_pot(x,y,z,rho_lm)
    enddo
  enddo
END SUBROUTINE

!----------------------------
FUNCTION mp_pot(x,y,z,rho_lm)
!----------------------------
  USE m_qd3d
  IMPLICIT NONE 
  integer :: L,lm
  REAL(8) :: mp_pot,xx,yy,zz,r,sump,x,y,z,rho_lm((L_max+1)**2)
  ! Function
  REAL(8) :: Ylm
  
  r = sqrt(x*x+y*y+z*z) + small
  xx = x/r
  yy = y/r
  zz = z/r
  sump = 0.d0
  do L = 0,L_max
    do lm=L**2+1,(L+1)**2
      sump = sump + Ylm(xx,yy,zz,lm)/r**(L+1)*rho_lm(lm)
    enddo
  enddo
  mp_pot = E2*sump
END FUNCTION

!---------------------
FUNCTION Ylm(x,y,z,lm)
!---------------------
  IMPLICIT NONE
	REAL(8), intent(IN) :: x,y,z
	integer, intent(IN) :: lm
	REAL(8)  :: Ylm,r
	! Spherical harmonics
	! It should be multiplied by r**l*sqrt((2*l+1)/4*pi)
	r=x**2+y**2+z**2
	select case( lm )
		case(1)  ; Ylm=1.d0                                     ! lm=1  (0  0)
		case(2)  ; Ylm=y                                        ! lm=2  (1 -1)
		case(3)  ; Ylm=z                                        ! lm=3  (1  0)
		case(4)  ; Ylm=x                                        ! lm=4  (1  1)
		case(5)  ; Ylm=sqrt(3.d0)*x*y                           ! lm=5  (2 -2)
		case(6)  ; Ylm=sqrt(3.d0)*y*z                           ! lm=6  (2 -1)
		case(7)  ; Ylm=(2*z*z-x*x-y*y)/2.d0                     ! lm=7  (2  0)
		case(8)  ; Ylm=sqrt(3.d0)*x*z                           ! lm=8  (2  1)
		case(9)  ; Ylm=sqrt(3.d0/4.d0)*(x*x-y*y)                ! lm=9  (2  2)
		case(10) ; Ylm=sqrt(5.d0/8.d0)*y*(3*x*x-y*y)            ! lm=10 (3 -3)
		case(11) ; Ylm=sqrt(15.d0)*x*y*z                        ! lm=11 (3 -2)
		case(12) ; Ylm=sqrt(3.d0/8.d0)*y*(4*z*z-x*x-y*y)        ! lm=12 (3 -1)
		case(13) ; Ylm=z*(2*z*z-3*x*x-3*y*y)/2.d0               ! lm=13 (3  0)
		case(14) ; Ylm=sqrt(3.d0/8.d0)*x*(4*z*z-x*x-y*y)        ! lm=14 (3  1)
		case(15) ; Ylm=sqrt(15.d0/4.d0)*z*(x*x-y*y)             ! lm=15 (3  2)
		case(16) ; Ylm=sqrt(5.d0/8.d0)*x*(x*x-3*y*y)            ! lm=16 (3  3)

		case(17) ; Ylm=sqrt(35.d0)/2.d0*x*y*(x**2-y**2)         ! lm=17 (4 -4)
		case(18) ; Ylm=sqrt(35.d0/8.d0)*y*z*(3*x**2-y**2)       ! lm=18 (4 -3)
		case(19) ; Ylm=sqrt(5.d0)/2.d0*x*y*(7*z**2-r**2)        ! lm=19 (4 -2)
		case(20) ; Ylm=sqrt(5.d0/8.d0)*y*(7*z**3-3*z*r**2)      ! lm=20 (4 -1)
		case(21) ; Ylm=(35*z**4-30*z**2*r**2+3.d0*r**2)/8.d0    ! lm=21 (4  0)
		case(22) ; Ylm=sqrt(5.d0/8.d0)*x*(7*z**3-3*z*r**2)      ! lm=22 (4  1)
		case(23) ; Ylm=sqrt(5.d0)/4.d0*(7*z**2-r**2)*(x**2-y**2)! lm=23 (4  2)
		case(24) ; Ylm=sqrt(35.d0/8.d0)*z*x*(x**2-3*y**2)       ! lm=24 (4  3)
		case(25) ; Ylm=sqrt(35.d0)/8.d0*(x**4+y**4-6*x**2*y**2) ! lm=25 (4  4)
	end select
END FUNCTION Ylm

!-----------------------------------
SUBROUTINE Hamiltonian_wavefn(ispin)
!-----------------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Calculate H|Phi>
  integer :: ispin
  if(N_d==4) then
	  call laplace_operator4()
	else
	  call laplace_operator()
	endif
  if(ispin==1) H_Phi = (V_ext+VH+V_exchange_up)*Phi - h2m*L_Phi
  if(ispin==2) H_Phi = (V_ext+VH+V_exchange_dw)*Phi - h2m*L_Phi
END SUBROUTINE 

!-----------------------------------
SUBROUTINE conjugate_gradient(ispin)
!-----------------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Conjugate gradient minimization to diagonalize the Hamiltonian
  integer :: orbital,iteration,i,ispin
  REAL(8)  :: Phi_H_Phi,gamma,beta_Phi,beta_beta,beta_H_Phi
  REAL(8)  :: Phi0_Phi0,Phi0_H_Phi0,delta,A,B,C,omega,overlap
  REAL(8), allocatable :: alpha(:),beta(:),Phi0(:),H_Phi0(:)

  allocate(alpha(N_L_points),beta(N_L_points),Phi0(N_L_points),H_Phi0(N_L_points))
  do orbital=1,N_orbitals(ispin)
		Phi=Psi(:,orbital,ispin)
		Phi0=Phi
		call Hamiltonian_wavefn(ispin)
		H_Phi0=H_Phi
		Phi0_H_Phi0=sum(Phi*H_Phi)*Grid_Volume
		Phi0_Phi0=1
		delta=Phi0_H_Phi0

    do iteration=1,N_iteration
			alpha=2*(H_Phi0-delta*Phi0)
			do i=1,orbital-1
				alpha=alpha-Psi(:,i,ispin)*sum(Psi(:,i,ispin)*alpha)*Grid_Volume
			enddo
			gamma=sum(alpha*alpha)*Grid_Volume
			if(iteration==1) beta=-alpha
			if(iteration>1) beta=-alpha+gamma/overlap*beta
			overlap=gamma
			beta_Phi=sum(Phi0*beta)*Grid_Volume
			beta_beta=sum(beta*beta)*Grid_Volume
			beta_H_Phi=sum(H_Phi0*beta)*Grid_Volume

			Phi=beta
			call Hamiltonian_wavefn(orbital)
			alpha=H_Phi
			Phi_H_Phi=sum(Phi*H_Phi)*Grid_Volume
			A = Phi_H_Phi*beta_Phi  - beta_H_Phi*beta_beta
			B = Phi_H_Phi*Phi0_Phi0   - Phi0_H_Phi0*beta_beta
			C = beta_H_Phi*Phi0_Phi0 - Phi0_H_Phi0*beta_Phi
			omega=(-B+sqrt(B*B-4*A*C))/(2*A)
			Phi0   = Phi0   +omega*Phi
			H_Phi0 = H_Phi0 +omega*H_Phi
			Phi0_Phi0   =sum(Phi0*Phi0)*Grid_Volume
			Phi0_H_Phi0 =sum(Phi0*H_Phi0)*Grid_Volume
			delta=Phi0_H_Phi0/Phi0_Phi0
		enddo
		Phi0_Phi0=sum(Phi0*Phi0)*Grid_Volume
		Psi(:,orbital,ispin)=Phi0/sqrt(Phi0_Phi0)
  enddo
  deallocate(alpha,beta,Phi0,H_Phi0)
END SUBROUTINE

!----------------------------------
SUBROUTINE orthogonalization(ispin)
!----------------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Gram-Schmidt orthogonalization
  integer :: orbital,i,ispin
  REAL(8) :: s
  complex(8) :: overlap ! FIXME: Need to be complex?

  do orbital=1,N_orbitals(ispin)
    do i=1,orbital-1
      overlap=sum(Psi(:,i,ispin)*Psi(:,orbital,ispin))*Grid_Volume
      Psi(:,orbital,ispin)=Psi(:,orbital,ispin)-overlap*Psi(:,i,ispin)
    enddo
    s=sum(abs(Psi(:,orbital,ispin))**2)*Grid_Volume
    Psi(:,orbital,ispin)=Psi(:,orbital,ispin)/sqrt(s)
  enddo
END SUBROUTINE 


!-------------------------------
SUBROUTINE Hamiltonian_density()
!-------------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Calculate the density dependent part of the Hamiltonian
  rho=density
  call bc_multipole()
  call poisson_solve_3d()
  vh = -v_pot
  call Exchange_Correlation()
END SUBROUTINE 


!---------------------------------
SUBROUTINE total_energy(iteration)
!---------------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Calculate the total energy
  integer :: i,orbital,ispin,iteration
  REAL(8)  :: E_sp,E_total                 
  
  E_sp=0.d0
  call Hamiltonian_density
  do ispin=1,2
    do orbital=1,N_orbitals(ispin)   
      Phi=Psi(:,orbital,ispin)
      call Hamiltonian_wavefn(ispin)
      Sp_Energy(orbital)=sum(Phi*H_Phi)*Grid_Volume
      E_sp=E_sp+Sp_energy(orbital)
    enddo
  enddo
  call Exchange_Correlation
  call Hartree_Energy
  write(*,*) 'E_sp    ', E_sp
  write(*,*) 'E_ex    ', E_exchange
  write(*,*) 'E_H     ', E_hartree
  E_total = E_sp + E_Exchange + E_Hartree
  write(*,*) 'E_total ', E_total
  
  IF( iteration > 0 ) THEN
    WRITE(2,'(i10,4ES16.8)') iteration, E_sp, E_exchange, E_hartree, E_total
  ENDIF 
END SUBROUTINE 

!----------------------------------
SUBROUTINE calculate_density(c1,c2)
!----------------------------------
  USE m_qd3d
  IMPLICIT NONE
  ! Calculate the density using linear mixing
  integer :: orbital,i,k,grid,ispin
  REAL(8)  :: c1,c2
		    
  Density_old=Density; Density=0.d0
  Density_up_old=Density_up; Density_dw_old=Density_dw
  Density_up=0.d0; Density_dw=0.d0 

  do ispin=1,2
  	do orbital=1,N_orbitals(ispin)
  	Phi(:)=Psi(:,orbital,ispin)
  	  do i=1,N_L_points
  	    if(ispin==1) Density_dw(i)=Density_dw(i)+Phi(i)*Phi(i)
  	    if(ispin==2) Density_up(i)=Density_up(i)+Phi(i)*Phi(i)
  	  enddo
  	enddo
  enddo
  Density_up=c1*Density_up_old+c2*Density_up
  Density_dw=c1*Density_dw_old+c2*Density_dw
  Density=Density_up+Density_dw
END SUBROUTINE

!--------------------------------
SUBROUTINE Exchange_Correlation()
!--------------------------------
  USE m_qd3d
  USE XC, ONLY: xc_pot
  IMPLICIT NONE

  ! Calculate the XC potential and energy
  integer :: i
  REAL(8)  :: su,eps_c,eps_x,rho_up,rho_down,v_x_up,v_x_down,v_c_up,v_c_down

  su=0.d0
  do i=1,N_L_points
    rho_up=density_up(i)+1.d-20
    rho_down=density_dw(i)+1.d-20
    call xc_pot(rho_up,rho_down,v_x_up,v_x_down,v_c_up,v_c_down,eps_x,eps_c)
    V_Exchange_up(i)=2*Ry*(v_x_up+v_c_up)
    V_Exchange_dw(i)=2*Ry*(v_x_down+v_c_down)
    su=su+Density(i)*(eps_x+eps_c)-Density_up(i)*V_exchange_up(i)-Density_dw(i)*V_exchange_dw(i)
  enddo
  E_exchange=su*grid_volume
END SUBROUTINE

!--------------------------
SUBROUTINE Hartree_Energy()
!--------------------------
  USE m_qd3d
  ! Calculate the Hartree-Energy
  E_Hartree=-0.5d0*dot_product(density,VH)*grid_volume
END SUBROUTINE Hartree_Energy

!------------------------------------------------------------------------------
PROGRAM qd3d
!------------------------------------------------------------------------------
  USE m_qd3d
  implicit none
  integer :: i,k,i1,i2,i3,num_args,ispin,kk,orbital,slice1,slice2,slice3
  REAL(8) :: x,y,z,x0,y0,z0,r,energy,the_sum
  REAL(8), allocatable :: dx_up(:),dy_up(:),dz_up(:),dx_dw(:),dy_dw(:),dz_dw(:)
  character(255) :: temp,total_density_filename,energy_vs_iteration_filename
  character(255) :: confining_potential_filename
  
  total_density_filename = "total_3d_density.3D.dat"
  confining_potential_filename = "confining_potential.3D.dat"
  energy_vs_iteration_filename = "energy_vs_iteration.dat"
  
  ! INPUT-----------------------------------
  N_orbitals(1) = 1
  N_orbitals(2) = 1
  omega0 = 2.d0
  N_L = (/ 20, 20, 20 /)
  grid_step = (/ 0.3d0, 0.3d0, 0.3d0 /)
  !-----------------------------------------

  N_L_points = product(N_L)
  grid_volume = product(grid_step)
  allocate(sp_energy((N_orbitals(1)+N_orbitals(2))),Psi(N_L_points,(N_orbitals(1)+N_orbitals(2)),2))
  allocate(Lattice(3,N_L_points),Lattice_inv(N_L(1),N_L(2),N_L(3)),grid_point(3,N_L_Points))
  allocate(V_X0(N_L(2),N_L(3)),V_XN(N_L(2),N_L(3)),V_Y0(N_L(1),N_L(3)),V_YN(N_L(1),N_L(3)))
  allocate(V_Z0(N_L(1),N_L(2)),V_ZN(N_L(1),N_L(2)),rho(N_L_points),V_POT(N_L_points))
  allocate(density(N_L_points),density_old(N_L_points),density_up(N_L_points))
  allocate(density_up_old(N_L_points),density_dw(N_L_points),density_dw_old(N_L_points))
  allocate(phi(N_L_points),L_phi(N_L_points),VH(N_L_points),V_ext(N_L_points))
  allocate(V_exchange_up(N_L_points),V_exchange_dw(N_L_points),H_Phi(N_L_Points))
  allocate(wf(-N_d:N_L(1)+N_d,-N_d:N_L(2)+N_d,-N_d:N_L(3)+N_d))  
  allocate(dx_up(N_L(1)),dy_up(N_L(2)),dz_up(N_L(3)),dx_dw(N_L(1)),dy_dw(N_L(2)),dz_dw(N_L(3)))
    
  ! Setup the lattice and initial guess for wavefunctions  
  call init_lattice
  do ispin=1,2
    do k=1,N_orbitals(ispin)
      call random_number(x0)
      call random_number(y0)
      call random_number(z0)
      do i=1,N_L_points
        x=grid_point(1,i)-x0
        y=grid_point(2,i)-y0
        z=grid_point(3,i)-z0
        Psi(i,k,ispin)=exp(-0.5d0*omega0*(x**2+y**2+z**2))
      enddo
    enddo
  enddo
  
  call init_confining_potential()
  
  ! Orthogonalize the initial wavefunctions, and use them to calculate the initial density and energy
  call orthogonalization(1)
  call orthogonalization(2)
  call calculate_density(0.d0,1.d0)
  call total_energy(0)

  ! Use the conjugate gradient method to diagonalize the Hamiltonian
  open(2,file=trim(adjustl(energy_vs_iteration_filename)))
  write(2,*)"# iteration E_sp E_exchange E_hartree E_total"
  do k=1,N_scf_iter
    write(*,*) k
    call orthogonalization(1)
    call conjugate_gradient(1)
    call orthogonalization(2)
    call conjugate_gradient(2)
    call calculate_density(0.5d0,0.5d0)
    call total_energy(k)
  enddo  
  close(2)
    
  ! Output the final 3D total density to a Point3D format file, suitable for VisIt
  open(1,file=trim(adjustl(total_density_filename)))
  open(2,file=trim(adjustl(confining_potential_filename)))
  write(1,*)"x y z density"
  write(2,*)"x y z V_ext"
  do i=1,N_L_points
    write(1,'(4ES16.8)')grid_point(1,i),grid_point(2,i),grid_point(3,i),density(i)*grid_volume
    write(2,'(4ES16.8)')grid_point(1,i),grid_point(2,i),grid_point(3,i),V_ext(i)   
  enddo
  close(1)
  close(2)
      
  deallocate(V_X0,V_XN,V_Y0,V_YN,V_Z0,V_ZN,rho,V_POT,density,density_old,density_up)
  deallocate(V_exchange_up,V_exchange_dw,H_Phi,dx_up,dy_up,dz_up,dx_dw,dy_dw,dz_dw)
  deallocate(density_up_old,density_dw,density_dw_old,phi,L_phi,VH,V_ext,wf)
  deallocate(sp_energy,Psi,lattice,lattice_inv,grid_point)
END PROGRAM 

