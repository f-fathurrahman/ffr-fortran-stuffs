!=========================================================================
subroutine write_energy_qp(nstate,energy_qp)
 use m_definitions
 use m_mpi
 use m_inputparam,only: nspin
 implicit none

 integer,intent(in)  :: nstate
 real(dp),intent(in) :: energy_qp(nstate,nspin)
!=====
 integer           :: energy_qpfile
 integer           :: istate
!=====

 !
 ! Only the proc iomaster writes down the ENERGY_QP file
 if( .NOT. is_iomaster) return

 write(stdout,'(/,a)') ' Writing ENERGY_QP file'


 open(newunit=energy_qpfile,file='ENERGY_QP',form='formatted')

 write(energy_qpfile,*) nspin
 write(energy_qpfile,*) nstate
 do istate=1,nstate
   write(energy_qpfile,*) istate,energy_qp(istate,:)
 enddo

 close(energy_qpfile)


end subroutine write_energy_qp