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