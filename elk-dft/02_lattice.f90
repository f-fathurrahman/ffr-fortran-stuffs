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

!-----------------------
subroutine prepare_all()
!-----------------------
  use modmpi
  use modomp
  use modmain, only: task, tasks
  use modvars, only: writevars
  implicit none

  call elk_init()

  ! read input files
  call readinput()

  ! initialise OpenMP variables
  call omp_init()

  ! initialise the MKL library
  call mkl_init()

  ! initialise the OpenBLAS library
  call oblas_init()

  ! initialise the BLIS library
  call blis_init()

  ! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)

  task = tasks(1)
  ! write task to VARIABLES.OUT
  call writevars('task',iv=task)

  call init0()
  call init1()

end subroutine

!------------------------------------------------------------------------------
program test
!------------------------------------------------------------------------------
  use modmain
  implicit none
  real(8) :: res(3,3)

  call prepare_all()

  write(*,*)
  write(*,*) 'avec'
  write(*,*)
  write(*,'(1x,3F18.10)') avec(1,:)
  write(*,'(1x,3F18.10)') avec(2,:)
  write(*,'(1x,3F18.10)') avec(3,:)

  write(*,*)
  write(*,*) 'bvec'
  write(*,*)
  write(*,'(1x,3F18.10)') bvec(1,:)
  write(*,'(1x,3F18.10)') bvec(2,:)
  write(*,'(1x,3F18.10)') bvec(3,:)

  res = matmul(avec,transpose(bvec))/(twopi)

  write(*,*)
  write(*,*) 'res = avec*bvec^T'
  write(*,*)
  write(*,'(1x,3F18.10)') res(1,:)
  write(*,'(1x,3F18.10)') res(2,:)
  write(*,'(1x,3F18.10)') res(3,:)


  call elk_stop()
end program
