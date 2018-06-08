!=========================================================================
subroutine set_occupation(nstate,temperature,electrons,magnetization,energy,occupation)
 use m_definitions
 use m_inputparam,only: print_matrix_, spin_fact, nspin
 use m_warning
 implicit none
 integer,intent(in)   :: nstate
 real(dp),intent(in)  :: electrons,magnetization,temperature
 real(dp),intent(in)  :: energy(nstate,nspin)
 real(dp),intent(out) :: occupation(nstate,nspin)
!=====
 real(dp)             :: remaining_electrons(nspin)
 real(dp)             :: electrons_mu,mu,delta_mu,grad_electrons
 real(dp)             :: mu_change
 integer              :: istate,nlines,ilines
 logical              :: file_exists
 integer              :: occfile
 integer              :: iter
!=====

 if( temperature < 1.0e-8_dp ) then

   occupation(:,:)=0.0_dp
  
   inquire(file='manual_occupations',exist=file_exists)
  
   if(.NOT. file_exists) then
     remaining_electrons(1) = (electrons+magnetization) / REAL(nspin,dp)
     if(nspin==2) remaining_electrons(2) = (electrons-magnetization) / REAL(nspin,dp)
  
     do istate=1,nstate
       occupation(istate,:) = MIN(remaining_electrons(:), spin_fact)
       remaining_electrons(:)  = remaining_electrons(:) - occupation(istate,:)
     enddo
   else
     write(stdout,*)
     write(stdout,*) 'occupations are read from file: manual_occupations'
     msg='reading manual occupations from file'
     call issue_warning(msg)
     open(newunit=occfile,file='manual_occupations',status='old')
     !
     ! read nlines, all other occupations are set to zero
     read(occfile,*) nlines
     do ilines=1,nlines
       read(occfile,*) occupation(ilines,:)
     enddo
     close(occfile)
     write(stdout,*) 'occupations set, closing file'
   endif

 else

   !
   ! Finite temperature case
   !
   write(stdout,'(1x,a,f12.6)') 'Find new the occupations and Fermi level for temperature (Ha): ',temperature

   mu = -0.1_dp
   delta_mu = 1.0e-5_dp
   electrons_mu = -1.0_dp
   iter = 0
   mu_change = 0.0_dp

   do while( ABS( electrons_mu - electrons ) > 1.0e-8_dp .AND. iter <= 100 )

     iter = iter + 1
     mu = mu + mu_change

     occupation(:,:) = fermi_dirac(energy,mu)
     electrons_mu    = SUM( occupation(:,:) )

     grad_electrons = ( SUM( fermi_dirac(energy,mu+delta_mu) ) - SUM( fermi_dirac(energy,mu-delta_mu) ) ) / ( 2.0_dp* delta_mu )

     ! Maximum change is made bounded within +/- 0.10 Hartree
     mu_change = -( electrons_mu - electrons ) / grad_electrons
     mu_change = MAX( MIN( mu_change , 0.10_dp / REAL(iter) ), -0.10_dp / REAL(iter) )

!     write(*,*) iter,mu,mu_change,0.10_dp / REAL(iter),electrons_mu

   enddo

   write(stdout,'(1x,a,f12.6)') 'Fermi level (eV): ', mu * Ha_eV

 endif

 !
 ! final check
 if( ABS( SUM(occupation(:,:)) - electrons ) > 1.0e-7_dp ) then
   write(stdout,*) 'occupation set up failed to give the right number of electrons'
   write(stdout,*) 'sum of occupations',SUM(occupation(:,:))
   write(stdout,*) 'electrons',electrons
   do istate=1,nstate
     write(stdout,*) istate,occupation(istate,:)
   enddo
   call die('FAILURE in set_occupation')
 endif 

 call dump_out_occupation('=== Occupations ===',nstate,nspin,occupation)

contains

function fermi_dirac(energy,mu)
 implicit none
 real(dp),intent(in) :: energy(nstate,nspin)
 real(dp),intent(in) :: mu
 real(dp)            :: fermi_dirac(nstate,nspin)
!=====

 fermi_dirac(:,:) = spin_fact / ( 1.0_dp + EXP( ( energy(:,:) - mu ) / temperature ) )

end function fermi_dirac

end subroutine set_occupation
