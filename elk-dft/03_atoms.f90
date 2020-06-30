include 'elk_init.f90'
include 'elk_stop.f90'
include 'prepare_all.f90'

!------------------------------------------------------------------------------
program test
!------------------------------------------------------------------------------
  use modmain
  implicit none
  integer :: is, ia

  call prepare_all()

  write(*,*) 'nspecies = ', nspecies
  write(*,*) 'natoms   = ', natoms(1:nspecies)
  write(*,*) 'natmmax  = ', natmmax

  write(*,*)
  write(*,*) 'atomic in lattice coordinates'
  write(*,*)
  do is=1,nspecies
    write(*,*) 'Species ', is, trim(spfname(is))
    do ia=1,natoms(is)
      write(*,*) 'atposl   = ', atposl(1:3,ia,is)
    enddo
  enddo

  write(*,*)
  write(*,*) 'atomic in Cartesian coordinates'
  write(*,*)
  do is=1,nspecies
    write(*,*) 'Species ', is, trim(spfname(is))
    do ia=1,natoms(is)
      write(*,*) 'atposc   = ', atposc(1:3,ia,is)
    enddo
  enddo

  call elk_stop()
end program
