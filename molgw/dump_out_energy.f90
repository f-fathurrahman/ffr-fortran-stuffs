!=========================================================================
subroutine dump_out_energy(title,nstate,nspin,occupation,energy)
 use m_definitions
 use m_mpi
 implicit none
 character(len=*),intent(in) :: title
 integer,intent(in)          :: nstate,nspin
 real(dp),intent(in)         :: occupation(nstate,nspin),energy(nstate,nspin)
!=====
 integer,parameter :: MAXSIZE=300
!=====
 integer  :: istate,ispin
 real(dp) :: spin_fact
!=====

 spin_fact = REAL(-nspin+3,dp)

 write(stdout,'(/,1x,a)') TRIM(title)

 if(nspin==1) then
   write(stdout,'(a)') '   #       (Ha)         (eV)      '
 else
   write(stdout,'(a)') '   #              (Ha)                      (eV)      '
   write(stdout,'(a)') '           spin 1       spin 2       spin 1       spin 2'
 endif
 do istate=1,MIN(nstate,MAXSIZE)
   select case(nspin)
   case(1)
     write(stdout,'(1x,i3,2(1x,f12.5),4x,f8.4)') istate,energy(istate,:),energy(istate,:)*Ha_eV,occupation(istate,:)
   case(2)
     write(stdout,'(1x,i3,2(2(1x,f12.5)),4x,2(f8.4,2x))') istate,energy(istate,:),energy(istate,:)*Ha_eV,occupation(istate,:)
   end select
   if(istate < nstate) then
     if( ANY( occupation(istate+1,:) < spin_fact/2.0_dp .AND. occupation(istate,:) > spin_fact/2.0 ) ) then 
        if(nspin==1) then
          write(stdout,'(a)') '  -----------------------------'
        else
          write(stdout,'(a)') '  -------------------------------------------------------'
        endif
     endif
   endif
 enddo

 write(stdout,*)

end subroutine dump_out_energy