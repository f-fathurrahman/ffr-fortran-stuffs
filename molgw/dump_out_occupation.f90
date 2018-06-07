!=========================================================================
subroutine dump_out_occupation(title,nstate,nspin,occupation)
 use m_definitions
 use m_mpi
 implicit none
 character(len=*),intent(in) :: title
 integer,intent(in)          :: nstate,nspin
 real(dp),intent(in)         :: occupation(nstate,nspin)
!=====
 integer :: ihomo
 integer :: istate,ispin
!=====

 write(stdout,'(/,1x,a)') TRIM(title)

 if( nspin == 2 ) then
   write(stdout,'(a)') '           spin 1       spin 2 '
 endif
 do istate=1,nstate
   if( ANY(occupation(istate,:) > 0.001_dp) ) ihomo = istate 
 enddo

 do istate=MAX(1,ihomo-5),MIN(ihomo+5,nstate)
   write(stdout,'(1x,i3,2(2(1x,f12.5)),2x)') istate,occupation(istate,:)
 enddo
 write(stdout,*)

end subroutine dump_out_occupation
