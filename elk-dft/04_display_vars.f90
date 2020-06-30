include 'elk_init.f90'
include 'elk_stop.f90'
include 'prepare_all.f90'
include 'write_apwlo_vars.f90'

!------------------------------------------------------------------------------
program test
!------------------------------------------------------------------------------
  use modmain
  implicit none
  integer :: is, ia

  call prepare_all()

  call write_apwlo_vars()

  call elk_stop()
end program
