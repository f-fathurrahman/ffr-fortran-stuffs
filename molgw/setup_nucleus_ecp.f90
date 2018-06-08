!=========================================================================
subroutine setup_nucleus_ecp(print_matrix_,basis,hamiltonian_nucleus)
  use m_definitions
 use m_atoms
 use m_dft_grid
 use m_ecp
 use m_timing
 use m_mpi
 use m_basis_set
 implicit none
 logical,intent(in)         :: print_matrix_
 type(basis_set),intent(in) :: basis
 real(dp),intent(inout)     :: hamiltonian_nucleus(basis%nbf,basis%nbf)
!=====
 integer              :: natom_local
 integer              :: ibf,jbf
 integer              :: iatom
 real(dp)             :: vnucleus_ij

 integer              :: iecp
 integer              :: iproj,nproj
 integer              :: mm
 integer              :: igrid
 real(dp)             :: rr(3)
 real(dp)             :: weight
 real(dp)             :: basis_function_r(basis%nbf)
 integer              :: iradial
 integer              :: i1,n1
 real(dp)             :: xtmp,phi,cos_theta
 real(dp)             :: wxa(nradial_ecp),xa(nradial_ecp)
 real(dp)             :: w1(nangular_ecp),x1(nangular_ecp),y1(nangular_ecp),z1(nangular_ecp)
 real(dp),allocatable :: int_fixed_r(:,:)
 real(dp),external    :: real_spherical_harmonics
 integer              :: necp,ie
 character(len=100)   :: title
 logical              :: element_has_ecp
!=====
  integer :: number_basis_function_am

 ! Check if there are some ECP
 if( nelement_ecp == 0 ) return


 call start_clock(timing_ecp)

 !
 ! Since there will be an allreduce operation in the end, 
 ! anticipate by dividing the input value of Hnucl by the number of procs
 if( nproc_world > 1 ) then
   hamiltonian_nucleus(:,:) = hamiltonian_nucleus(:,:) / nproc_world
 endif

 n1 = nangular_ecp
 select case(nangular_ecp)
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
   write(stdout,*) 'grid points: ',nangular_ecp
   call die('setup_nucleus_ecp: Lebedev grid is not available')
 end select


 do iradial=1,nradial_ecp
   xtmp = ( iradial - 0.5_dp ) / REAL(nradial_ecp,dp)
   xa(iradial)   = -5.0_dp * log( 1.0_dp - xtmp**3)
   wxa(iradial)  = 3.0_dp * 5.0_dp * xtmp**2 / ( 1.0_dp - xtmp**3 ) / REAL(nradial_ecp,dp)
 enddo


 do iatom=1,natom
   element_has_ecp = .FALSE.
   do ie=1,nelement_ecp
     if( element_ecp(ie) == basis_element(iatom) ) then
       element_has_ecp = .TRUE.
       exit
     endif
   enddo 

   if( .NOT. element_has_ecp ) cycle

   necp = ecp(ie)%necp
     

   nproj = 0
   do iecp=1,necp
     nproj = nproj + number_basis_function_am('PURE',ecp(ie)%lk(iecp))
   enddo
  
  
   allocate(int_fixed_r(basis%nbf,nproj))
   do iradial=1,nradial_ecp
     if( MODULO(iradial-1,nproc_world) /= rank_world ) cycle

     int_fixed_r(:,:) = 0.0_dp
     do i1=1,nangular_ecp
       rr(1) = xa(iradial) * x1(i1) + x(1,iatom)
       rr(2) = xa(iradial) * y1(i1) + x(2,iatom)
       rr(3) = xa(iradial) * z1(i1) + x(3,iatom)
       call calculate_basis_functions_r(basis,rr,basis_function_r)
  
       cos_theta = z1(i1)
       phi       = ATAN2(y1(i1),x1(i1))
  
       iproj = 0
       do iecp=1,necp
         do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
           iproj = iproj + 1
           int_fixed_r(:,iproj) = int_fixed_r(:,iproj) + basis_function_r(:) &
                                     * real_spherical_harmonics(ecp(ie)%lk(iecp),mm,cos_theta,phi) &
                                        * w1(i1) * 4.0_dp * pi  
         enddo
       enddo
     enddo
  
  
     iproj = 0
     do iecp=1,necp
       do mm=-ecp(ie)%lk(iecp),ecp(ie)%lk(iecp)
         iproj = iproj + 1
         do jbf=1,basis%nbf
           do ibf=1,basis%nbf
             hamiltonian_nucleus(ibf,jbf) = hamiltonian_nucleus(ibf,jbf)  &
                 + int_fixed_r(ibf,iproj) * int_fixed_r(jbf,iproj) * wxa(iradial) * xa(iradial)**2  &
                    * ecp(ie)%dk(iecp) * EXP( -ecp(ie)%zetak(iecp) * xa(iradial)**2 ) * xa(iradial)**(ecp(ie)%nk(iecp)-2)
           enddo
         enddo
       enddo
     enddo
  
   enddo
  
   deallocate(int_fixed_r)

 enddo 

 call xsum_world(hamiltonian_nucleus)

 title='=== ECP Nucleus potential contribution ==='
 call dump_out_matrix(print_matrix_,title,basis%nbf,1,hamiltonian_nucleus)

 call stop_clock(timing_ecp)

end subroutine setup_nucleus_ecp
