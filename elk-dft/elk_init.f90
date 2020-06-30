!--------------------
subroutine elk_init()
!--------------------
  use modmpi
  use modmain, only: version
  implicit none
  
  ! initialise MPI execution environment
  call mpi_init(ierror)
  
  ! duplicate mpi_comm_world
  call mpi_comm_dup(mpi_comm_world,mpicom,ierror)
  
  ! determine the number of MPI processes
  call mpi_comm_size(mpicom,np_mpi,ierror)
  
  ! determine the local MPI process number
  call mpi_comm_rank(mpicom,lp_mpi,ierror)
  
  ! determine if the local process is the master
  if (lp_mpi.eq.0) then
    mp_mpi=.true.
    write(*,*)
    write(*,'("Elk code version ",I1.1,".",I1.1,".",I2.2," started")') version
  else
    mp_mpi=.false.
  end if
end subroutine