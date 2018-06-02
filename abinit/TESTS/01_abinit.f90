PROGRAM test

  USE defs_basis, ONLY : abi_io_redirect

  USE m_xmpi, ONLY : xmpi_world
  USE m_xmpi, ONLY : xmpi_init, xmpi_comm_rank, xmpi_end

  IMPLICIT NONE 
  INTEGER :: me
  TYPE(args_t) :: args

  CALL abi_io_redirect(new_io_comm=xmpi_world)

  CALL xmpi_init()
  me = xmpi_comm_rank(xmpi_world)

  WRITE(*,*) 'me = ', me

  ! Parse command line arguments.
  args = args_parser()
  IF( args%EXIT /= 0 ) GOTO 100

  WRITE(*,*)
  WRITE(*,*) 'Pass here'

100 CALL xmpi_end()

END PROGRAM 

