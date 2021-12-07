include 'elk_init.f90'
include 'elk_stop.f90'
include 'prepare_all.f90'
include 'info_lattice.f90'
include 'info_muffin_tins.f90'
include 'info_apwlo.f90'

!------------------------------------------------------------------------------
program test
!------------------------------------------------------------------------------
  use modmain
  implicit none
  integer :: is, ia, ir

  call prepare_all()

  call info_lattice()
  !call info_muffin_tins()
  !call info_apwlo()
  
  call elk_stop()
end program
