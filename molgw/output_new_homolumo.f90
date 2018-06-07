!=========================================================================
subroutine output_new_homolumo(calculation_name,nstate,occupation,energy,istate_min,istate_max)
 use m_definitions
 use m_mpi
 use m_inputparam,only: nspin,spin_fact
 implicit none

 character(len=*),intent(in) :: calculation_name
 integer,intent(in)          :: nstate,istate_min,istate_max
 real(dp),intent(in)         :: occupation(nstate,nspin),energy(nstate,nspin)
!=====
 real(dp) :: ehomo_tmp,elumo_tmp
 real(dp) :: ehomo(nspin),elumo(nspin)
 integer  :: ispin,istate
!=====

 do ispin=1,nspin
   ehomo_tmp=-HUGE(1.0_dp)
   elumo_tmp= HUGE(1.0_dp)

   do istate=istate_min,istate_max

     if( occupation(istate,ispin)/spin_fact > completely_empty ) then
       ehomo_tmp = MAX( ehomo_tmp , energy(istate,ispin) )
     endif

     if( occupation(istate,ispin)/spin_fact < 1.0_dp - completely_empty ) then
       elumo_tmp = MIN( elumo_tmp , energy(istate,ispin) )
     endif

   enddo

   ehomo(ispin) = ehomo_tmp
   elumo(ispin) = elumo_tmp

 enddo


 write(stdout,*)
 if( ALL( ehomo(:) > -1.0e6 ) ) then
   write(stdout,'(1x,a,1x,a,2(3x,f12.6))') TRIM(calculation_name),'HOMO energy    (eV):',ehomo(:) * Ha_eV
 endif
 if( ALL( elumo(:) <  1.0e6 ) ) then
   write(stdout,'(1x,a,1x,a,2(3x,f12.6))') TRIM(calculation_name),'LUMO energy    (eV):',elumo(:) * Ha_eV
 endif
 if( ALL( ehomo(:) > -1.0e6 ) .AND. ALL( elumo(:) <  1.0e6 ) ) then
   write(stdout,'(1x,a,1x,a,2(3x,f12.6))') TRIM(calculation_name),'HOMO-LUMO gap  (eV):',( elumo(:)-ehomo(:) ) * Ha_eV
 endif
 write(stdout,*)


end subroutine output_new_homolumo