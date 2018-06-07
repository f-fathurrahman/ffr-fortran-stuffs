!=========================================================================
subroutine init_basis_set(basis_path,basis_name,ecp_basis_name,gaussian_type,basis)
  use m_definitions
  use m_atoms
  use m_basis_set
  use m_ecp
 implicit none
 character(len=4),intent(in)   :: gaussian_type
 character(len=100),intent(in) :: basis_path
 character(len=100),intent(in) :: basis_name(natom_basis)
 character(len=100),intent(in) :: ecp_basis_name(natom_basis)
 type(basis_set),intent(out)   :: basis
  !=====
  character(len=100)            :: basis_filename
  integer                       :: ibf,jbf,kbf,ng,ig
  integer                       :: ishell,ishell_file
  integer                       :: jbf_cart
  real(dp),allocatable          :: alpha(:),coeff(:)
  logical                       :: file_exists
  integer                       :: basisfile
  integer                       :: am_read,nshell_file
  logical,parameter             :: normalized=.TRUE.
  integer                       :: iatom
  integer                       :: index_in_shell
  integer                       :: nx,ny,nz,mm
  real(dp)                      :: x0(3)
  !=====
  integer :: number_basis_function_am

 basis%nbf           = 0
 basis%nbf_cart      = 0
 basis%nshell        = 0
 basis%gaussian_type = gaussian_type

 if(TRIM(basis_name(1))=='none') return

 !
 ! LOOP OVER ATOMS
 !
 do iatom=1,natom_basis

   if( nelement_ecp > 0 ) then 
     if( ANY( element_ecp(:) == basis_element(iatom) ) ) then
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(ecp_basis_name(iatom)))
     else
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
     endif
   else
     basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
   endif
  
   inquire(file=TRIM(basis_filename),exist=file_exists)
   if(.NOT.file_exists) then
     write(stdout,'(a,a)') ' Looking for file ',TRIM(basis_filename)
     call die('basis set file not found')
   endif
  
   !
   ! read first to get all the dimensions
   open(newunit=basisfile,file=TRIM(basis_filename),status='old')
   read(basisfile,*) nshell_file
   if(nshell_file<1) call die('ERROR in basis set file')
   do ishell_file=1,nshell_file
     read(basisfile,*) ng,am_read
     if(ng<1) call die('ERROR in basis set file')
     if(am_read==10) call die('Deprecated basis set file with shared exponent SP orbitals. Please split them')
     basis%nbf_cart = basis%nbf_cart + number_basis_function_am('CART'             ,am_read)
     basis%nbf      = basis%nbf      + number_basis_function_am(basis%gaussian_type,am_read)
     basis%nshell   = basis%nshell   + 1
     do ig=1,ng
       read(basisfile,*) 
     enddo
   enddo
   close(basisfile)
  
 enddo


 write(stdout,*)
 write(stdout,'(a50,i8)') 'Total number of basis functions:',basis%nbf
 if(basis%gaussian_type=='PURE') then
   write(stdout,'(a50,i8)') 'Total number of cart. functions:',basis%nbf_cart
 endif
 write(stdout,'(a50,i8)') 'Number of shells:',basis%nshell

 allocate(basis%bfc(basis%nbf_cart))
 allocate(basis%bff(basis%nbf))
 allocate(basis%shell(basis%nshell))

 jbf         = 0
 jbf_cart    = 0
 ishell      = 0
 do iatom=1,natom_basis

   if( nelement_ecp > 0 ) then 
     if( ANY( element_ecp(:) == basis_element(iatom) ) ) then
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(ecp_basis_name(iatom)))
     else
       basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
     endif
   else
     basis_filename=ADJUSTL(TRIM(basis_path)//'/'//TRIM(ADJUSTL(element_name(REAL(basis_element(iatom),dp))))//'_'//TRIM(basis_name(iatom)))
   endif
  
   open(newunit=basisfile,file=TRIM(basis_filename),status='old')
   read(basisfile,*) nshell_file
   do ishell_file=1,nshell_file
     read(basisfile,*) ng,am_read
     allocate(alpha(ng),coeff(ng))
  
     do ig=1,ng
       read(basisfile,*) alpha(ig),coeff(ig)
     enddo
  
     x0(:) = x(:,iatom)

     !
     ! Shell setup
     !
     ishell = ishell + 1

     basis%shell(ishell)%ishell = ishell
     basis%shell(ishell)%am     = am_read
     basis%shell(ishell)%iatom  = iatom
     basis%shell(ishell)%x0(:)  = x0(:)
     basis%shell(ishell)%ng     = ng
     allocate(basis%shell(ishell)%alpha(ng))
     allocate(basis%shell(ishell)%coeff(ng))
     basis%shell(ishell)%alpha(:)    = alpha(:)
     basis%shell(ishell)%istart      = jbf + 1
     basis%shell(ishell)%iend        = jbf + number_basis_function_am(gaussian_type,am_read)
     basis%shell(ishell)%istart_cart = jbf_cart + 1
     basis%shell(ishell)%iend_cart   = jbf_cart + number_basis_function_am('CART',am_read)
     ! shell%coeff(:) is setup just after the basis functions


     !
     ! Basis function setup
     !

     !
     ! Ordering of Libint as explained in Kenny et al. J. Comput Chem. 29, 562 (2008).
     !
     nx = am_read
     ny = 0
     nz = 0
     index_in_shell = 0
     do 
       ! Add the new basis function
       jbf_cart = jbf_cart + 1 
       index_in_shell = index_in_shell + 1
       call init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,ishell,index_in_shell,basis%bfc(jbf_cart))
       if(basis%gaussian_type == 'CART') then
         jbf = jbf + 1
         call init_basis_function(normalized,ng,nx,ny,nz,iatom,x0,alpha,coeff,ishell,index_in_shell,basis%bff(jbf))
       endif

       ! Break the loop when nz is equal to l
       if( nz == am_read ) exit

       if( nz < am_read - nx ) then
         ny = ny - 1
         nz = nz + 1
       else
         nx = nx - 1
         ny = am_read - nx
         nz = 0
       endif

     enddo

     index_in_shell = 0
     if(basis%gaussian_type == 'PURE') then
       do mm=-am_read,am_read
         jbf = jbf + 1
         index_in_shell = index_in_shell + 1
         call init_basis_function_pure(normalized,ng,am_read,mm,iatom,x0,alpha,coeff,ishell,index_in_shell,basis%bff(jbf))
       enddo
     endif

     !
     ! Include here the normalization part that does not depend on (nx,ny,nz)
     basis%shell(ishell)%coeff(:) = basis%bfc(jbf_cart-number_basis_function_am(gaussian_type,am_read)+1)%coeff(:) &
               * ( 2.0_dp / pi )**0.75_dp * 2.0_dp**am_read * alpha(:)**( 0.25_dp * ( 2.0_dp*am_read + 3.0_dp ) )
  
     deallocate(alpha,coeff)
   enddo
   close(basisfile)

 !
 ! END OF THE LOOP OVER ATOMS
 enddo
 


 ! Find the maximum angular momentum employed in the basis set
 basis%ammax = MAXVAL(basis%bfc(:)%am)

 write(stdout,'(a50,i8)') 'Maximum angular momentum in the basis set:',basis%ammax
 write(stdout,'(a50,a8)') '                                          ',orbital_momentum_name(basis%ammax)

 if(basis%ammax > MOLGW_LMAX ) then      
   write(stdout,*) 'Maximum angular momentum: ',basis%ammax
   write(stdout,*) 'while this compilation of LIBINT only supports: ',MOLGW_LMAX
   call die('init_basis_set: Angular momentum is too high')
 endif

 !
 ! finally output the basis set for debugging
 if( .FALSE. ) then
   do ibf=1,basis%nbf_cart
     write(stdout,*) ' Cartesian function number',ibf
     call print_basis_function(basis%bfc(ibf))
   enddo
 endif

 write(stdout,'(a,/)') ' Basis set is fit and ready'

end subroutine init_basis_set