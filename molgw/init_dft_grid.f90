!=========================================================================
subroutine init_dft_grid(basis,grid_level_in,needs_gradient,precalculate_wfn,batch_size)
  use m_elements
  use m_atoms
  use m_tools,only: coeffs_gausslegint
  use m_basis_set
  use m_dft_grid
  use m_inputparam, only : partition_scheme, grid_memory
  use m_timing, only : timing_grid_generation
  USE m_mpi
  implicit none

 type(basis_set),intent(in) :: basis
 integer,intent(in)         :: grid_level_in
 logical,intent(in)         :: needs_gradient,precalculate_wfn
 integer,intent(in)         :: batch_size
  !=====
  integer              :: ngrid_max_allowed
  integer              :: iradial,iatom,iangular,ir,ir1,igrid
  integer              :: n1,n2,nangular,ngridmax
  real(dp)             :: weight,radius
  real(dp),allocatable :: x1(:),x2(:)
  real(dp),allocatable :: y1(:),y2(:)
  real(dp),allocatable :: z1(:),z2(:)
  real(dp),allocatable :: w1(:),w2(:)
  real(dp),allocatable :: xa(:,:),wxa(:,:)
  real(dp)             :: p_becke(natom_basis),s_becke(natom_basis,natom_basis),fact_becke 
  real(dp)             :: mu,alpha,xtmp,mu_aa
  integer              :: jatom,katom
  real(dp),allocatable :: rr_grid_tmp(:,:)
  real(dp),allocatable :: w_grid_tmp(:)
  !=====
  real(dp) :: smooth_step

 call start_clock(timing_grid_generation)
 ngrid_stored = 0
 select case(grid_level_in)
 case(low)       ! accuracy not guaranted, just for quick test runs
   nradial         =  25
   nangular_fine   =  26
   nangular_coarse =   6
 case(medium)    ! 10 meV accuracy on potentials
   nradial         =  40
   nangular_fine   =  50
   nangular_coarse =  14
 case(high)      !  1 meV accuracy on potentials
   nradial         =  60
   nangular_fine   = 110
   nangular_coarse =  38
 case(very_high) ! almost perfect potentials
   nradial         =  70 
   nangular_fine   = 170 
   nangular_coarse =  50 
 case(insane)    ! overdoing a lot
   nradial         = 200
   nangular_fine   = 434
   nangular_coarse = 434
 case default
   call die('integration quality not recognized')
 end select

 allocate(xa(nradial,natom_basis),wxa(nradial,natom_basis))
 allocate(x1(nangular_fine),y1(nangular_fine),z1(nangular_fine),w1(nangular_fine))
 allocate(x2(nangular_coarse),y2(nangular_coarse),z2(nangular_coarse),w2(nangular_coarse))

 !
 ! spherical integration
 !

! !
! ! radial part with Gauss-Legendre
! call coeffs_gausslegint(-1.0_dp,1.0_dp,xa(:,1),wxa(:,1),nradial)
! !
! ! Transformation from [-1;1] to [0;+\infty[
! ! taken from M. Krack JCP 1998
! wxa(:,1) = wxa(:,1) * ( 1.0_dp / log(2.0_dp) / ( 1.0_dp - xa(:,1) ) )
! xa(:,1)  = log( 2.0_dp / (1.0_dp - xa(:,1) ) ) / log(2.0_dp)
!
! ! All atoms have the same grid
! do iatom=2,natom_basis
!   xa(:,iatom)  =  xa(:,1)
!   wxa(:,iatom) = wxa(:,1)
! enddo


 !
 ! LOG3 radial grid
 do iatom=1,natom_basis

   select case(NINT(zatom(iatom)))
   case(3,4,11,12,19,20)
     alpha = 7.0_dp
   case default
     alpha = 5.0_dp
   end select

   do iradial=1,nradial
     xtmp = ( iradial - 0.5_dp ) / REAL(nradial,dp)
     xa(iradial,iatom)   = -alpha * log( 1.0_dp - xtmp**3)
     wxa(iradial,iatom)  = 3.0_dp * alpha * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nradial,dp)
   enddo

 enddo



 n1 = nangular_fine
 ! angular part with Lebedev - Laikov
 ! (x1,y1,z1) on the unit sphere
 select case(nangular_fine)
 case(6)
   call ld0006(x1,y1,z1,w1,n1)
 case(14)
   call ld0014(x1,y1,z1,w1,n1)
 case(26)
   call ld0026(x1,y1,z1,w1,n1)
 case(38)
   call ld0038(x1,y1,z1,w1,n1)
 case(50)
   call ld0050(x1,y1,z1,w1,n1)
 case(74)
   call ld0074(x1,y1,z1,w1,n1)
 case(86)
   call ld0086(x1,y1,z1,w1,n1)
 case(110)
   call ld0110(x1,y1,z1,w1,n1)
 case(146)
   call ld0146(x1,y1,z1,w1,n1)
 case(170)
   call ld0170(x1,y1,z1,w1,n1)
 case(230)
   call ld0230(x1,y1,z1,w1,n1)
 case(302)
   call ld0302(x1,y1,z1,w1,n1)
 case(434)
   call ld0434(x1,y1,z1,w1,n1)
 case default
   write(stdout,*) 'grid points: ',nangular_fine
   call die('Lebedev grid is not available')
 end select

 n2 = nangular_coarse
 ! angular part with Lebedev - Laikov
 ! (x2,y2,z2) on the unit sphere
 select case(nangular_coarse)
 case(6)
   call ld0006(x2,y2,z2,w2,n2)
 case(14)
   call ld0014(x2,y2,z2,w2,n2)
 case(26)
   call ld0026(x2,y2,z2,w2,n2)
 case(38)
   call ld0038(x2,y2,z2,w2,n2)
 case(50)
   call ld0050(x2,y2,z2,w2,n2)
 case(74)
   call ld0074(x2,y2,z2,w2,n2)
 case(86)
   call ld0086(x2,y2,z2,w2,n2)
 case(110)
   call ld0110(x2,y2,z2,w2,n2)
 case(146)
   call ld0146(x2,y2,z2,w2,n2)
 case(170)
   call ld0170(x2,y2,z2,w2,n2)
 case(230)
   call ld0230(x2,y2,z2,w2,n2)
 case(434)
   call ld0434(x2,y2,z2,w2,n2)
 case default
   write(stdout,*) 'grid points: ',nangular_coarse
   call die('Lebedev grid is not available')
 end select

 ! Calculate the maximum number of grid points
 ngridmax = 0
 do iatom=1,natom_basis
   radius = element_covalent_radius(REAL(basis_element(iatom),dp))
   do iradial=1,nradial
     if( xa(iradial,iatom) < pruning_limit * radius ) then
       ngridmax = ngridmax + nangular_coarse
     else
       ngridmax = ngridmax + nangular_fine
     endif
   enddo
 enddo

 !
 ! Temporary storage before the screening of the low weights
 allocate(rr_grid_tmp(3,ngridmax),w_grid_tmp(ngridmax))

 rr_grid_tmp(:,:) = 0.0_dp 
 w_grid_tmp(:)    = 0.0_dp
 ir    = 0
 do iatom=1,natom_basis
   radius = element_covalent_radius(REAL(basis_element(iatom),dp))

   do iradial=1,nradial
     if( xa(iradial,iatom) < pruning_limit * radius ) then
       nangular = nangular_coarse
     else
       nangular = nangular_fine
     endif

     do iangular=1,nangular
       ir = ir + 1

       ! Parallelization of the weights generation
       if( MODULO(ir,nproc_world) /= rank_world ) cycle

       if( xa(iradial,iatom) < pruning_limit * radius ) then
         rr_grid_tmp(1,ir) = xa(iradial,iatom) * x2(iangular) + x(1,iatom)
         rr_grid_tmp(2,ir) = xa(iradial,iatom) * y2(iangular) + x(2,iatom)
         rr_grid_tmp(3,ir) = xa(iradial,iatom) * z2(iangular) + x(3,iatom)
         weight   = wxa(iradial,iatom) * w2(iangular) * xa(iradial,iatom)**2 * 4.0_dp * pi
       else
         rr_grid_tmp(1,ir) = xa(iradial,iatom) * x1(iangular) + x(1,iatom)
         rr_grid_tmp(2,ir) = xa(iradial,iatom) * y1(iangular) + x(2,iatom)
         rr_grid_tmp(3,ir) = xa(iradial,iatom) * z1(iangular) + x(3,iatom)
         weight   = wxa(iradial,iatom) * w1(iangular) * xa(iradial,iatom)**2 * 4.0_dp * pi
       endif


       select case(partition_scheme)
       case('becke')
         !
         ! Partitionning scheme of Axel Becke, J. Chem. Phys. 88, 2547 (1988).
         !
         s_becke(:,:) = 0.0_dp
         do katom=1,natom_basis
           do jatom=1,natom_basis
             if(katom==jatom) cycle
             mu = ( NORM2(rr_grid_tmp(:,ir)-x(:,katom)) - NORM2(rr_grid_tmp(:,ir)-x(:,jatom)) ) &
                       / NORM2(x(:,katom)-x(:,jatom))
             s_becke(katom,jatom) = 0.5_dp * ( 1.0_dp - smooth_step(smooth_step(smooth_step(mu))) )
           enddo
         enddo

       case('ssf')
         !
         ! Partitionning scheme of Stratmann, Scuseria, Frisch, Chem. Phys. Lett. 257, 213 (1996)
         !
         s_becke(:,:) = 0.0_dp
         do katom=1,natom_basis
           do jatom=1,natom_basis
             if(katom==jatom) cycle
             mu = ( NORM2(rr_grid_tmp(:,ir)-x(:,katom)) - NORM2(rr_grid_tmp(:,ir)-x(:,jatom)) ) &
                       / NORM2(x(:,katom)-x(:,jatom))

             if( mu < -aa ) then
               s_becke(katom,jatom) = 1.0_dp
             else if( mu > aa ) then
               s_becke(katom,jatom) = 0.0_dp
             else 
               mu_aa = mu / aa
               s_becke(katom,jatom) = ( 35.0_dp * mu_aa - 35.0_dp * mu_aa**3 + 21.0_dp * mu_aa**5 - 5.0_dp * mu_aa**7 ) / 16.0_dp

               s_becke(katom,jatom) = 0.5_dp * ( 1.0_dp - s_becke(katom,jatom) )
             endif

           enddo
         enddo

       case default
         call die('Invalid choice for the partition scheme. Change partion_scheme value in the input file')
       end select

       p_becke(:) = 1.0_dp
       do katom=1,natom_basis
         do jatom=1,natom_basis
           if(katom==jatom) cycle
           p_becke(katom) = p_becke(katom) * s_becke(katom,jatom)
         enddo
       enddo
       fact_becke = p_becke(iatom) / SUM( p_becke(:) )

       w_grid_tmp(ir) = weight * fact_becke

     enddo
   enddo
 enddo

 deallocate(x1,y1,z1,w1)
 deallocate(x2,y2,z2,w2)

 call xsum_world(rr_grid_tmp)
 call xsum_world(w_grid_tmp)


 ! Denombrate the number of non-negligible weight in the quadrature
 ngrid = COUNT( w_grid_tmp(:) > TOL_WEIGHT )

 ! Distribute the grid over processors
 call init_dft_grid_distribution(ngrid)

 write(stdout,'(/,a)')              ' Setup the DFT quadrature'
 write(stdout,'(a,i4,1x,i4,1x,i4)') ' discretization grid per atom [radial , angular max - angular min] ',nradial,nangular_fine,nangular_coarse
 write(stdout,'(a,i8)')             ' total number of real-space points for this processor',ngrid

 allocate(rr_grid(3,ngrid),w_grid(ngrid))

 igrid = 0
 ir = 0
 do ir1=1,ngridmax
   if( w_grid_tmp(ir1) > TOL_WEIGHT ) then
     igrid = igrid + 1
     if( is_my_grid_task(igrid) ) then
       ir = ir + 1
       w_grid(ir) = w_grid_tmp(ir1)
       rr_grid(:,ir) = rr_grid_tmp(:,ir1)
     endif

   endif
 enddo

 deallocate(rr_grid_tmp,w_grid_tmp)


 !
 ! Once the grid is set up, 
 ! precalculate the wavefunctions and their gradient on a part of the grid
 !

 ! Save an internal copy of batch_size
 batch_size_ = batch_size


 if( precalculate_wfn ) then
   !
   ! grid_memory is given in Megabytes
   ! If gradient is needed, the storage is 4 times larger
   if( .NOT. needs_gradient ) then
     ngrid_max_allowed = NINT( grid_memory * 1024.0_dp**2 / ( REAL(basis%nbf,dp) * REAL(dp,dp) ) )
   else
     ngrid_max_allowed = NINT( grid_memory * 1024.0_dp**2 / ( 4.0_dp * REAL(basis%nbf,dp) * REAL(dp,dp) ) )
   endif

   ngrid_stored = MIN(ngrid,ngrid_max_allowed)
   ! Enforce a multiple of batches
   ngrid_stored = batch_size * ( ngrid_stored/batch_size )

   call prepare_basis_functions_r(basis,batch_size_)
   if( needs_gradient ) then
     call prepare_basis_functions_gradr(basis,batch_size_)
   endif

 else
   ngrid_stored = 0
 endif


 call stop_clock(timing_grid_generation)


end subroutine init_dft_grid
