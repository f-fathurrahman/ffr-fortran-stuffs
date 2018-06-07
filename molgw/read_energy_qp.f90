!=========================================================================
subroutine read_energy_qp(nstate,energy_qp,reading_status)
 use m_definitions
 use m_mpi
 use m_warning,only: issue_warning
 use m_inputparam,only: nspin
 implicit none

 integer,intent(in)   :: nstate
 integer,intent(out)  :: reading_status
 real(dp),intent(out) :: energy_qp(nstate,nspin)
!=====
 integer           :: energy_qpfile
 integer           :: istate,jstate
 integer           :: nspin_read,nstate_read
 logical           :: file_exists_capitalized,file_exists
!=====

 write(stdout,'(/,a)') ' Reading ENERGY_QP file'

 inquire(file='ENERGY_QP',exist=file_exists_capitalized)
 inquire(file='energy_qp',exist=file_exists)

 if(file_exists_capitalized) then
   open(newunit=energy_qpfile,file='ENERGY_QP',form='formatted',status='old')
 else if(file_exists) then
   open(newunit=energy_qpfile,file='energy_qp',form='formatted',status='old')
 endif

 if( file_exists_capitalized .OR. file_exists ) then
   read(energy_qpfile,*) nspin_read
   read(energy_qpfile,*) nstate_read
   if( nstate_read /= nstate .OR. nspin_read /= nspin ) then
     call issue_warning('ENERGY_QP file does not have the correct dimensions')
     reading_status=2
   else
     do istate=1,nstate
       read(energy_qpfile,*) jstate,energy_qp(istate,:)
       ! Scissor operator
       if( jstate == -1 ) then
         reading_status=-1
         close(energy_qpfile)
         return
       endif
     enddo
     reading_status=0
   endif
   close(energy_qpfile)
 else
   reading_status=1
   call issue_warning('files ENERGY_QP and energy_qp do not exist')
 endif


end subroutine read_energy_qp