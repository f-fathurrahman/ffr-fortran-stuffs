!--------------------
subroutine elk_stop()
!--------------------
  use modmpi
  implicit none

  if (mp_mpi) then
    open(50,file='RUNNING')
    close(50,status='DELETE')
    write(*,*)
    write(*,'("Elk code stopped")')
  end if
  ! terminate MPI execution environment
  call mpi_finalize(ierror)
end subroutine