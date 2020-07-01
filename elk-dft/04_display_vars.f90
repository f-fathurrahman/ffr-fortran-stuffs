include 'elk_init.f90'
include 'elk_stop.f90'
include 'prepare_all.f90'
include 'write_apwlo_vars.f90'

!------------------------------------------------------------------------------
program test
!------------------------------------------------------------------------------
  use modmain
  implicit none
  integer :: is, ia, ir

  call prepare_all()

  !call write_apwlo_vars()
  
  do is = 1,nspecies
    write(*,*) 'Species = ', trim(spfname(is))
    write(*,*) 'nrsp = ', nrsp(is)
    write(*,*) 'nrmt = ', nrmt(is)
    write(*,*) 'rminsp = ', rminsp(is)
    write(*,*) 'rmaxsp = ', rmaxsp(is)
    
    write(*,*)
    write(*,'(1x,A,5ES18.10)') 'wrmt = ', wrmt(1:5,is)
    write(*,'(1x,A,5ES18.10)') 'wrcmt = ', wrcmt(1:5,is)

    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprmt(1,1:5,is)
    !write(*,*)
    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprmt(2,1:5,is)
    !write(*,*)
    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprmt(3,1:5,is)
    !write(*,*)
    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprmt(4,1:5,is)

    write(*,'(1x,A,5ES18.10)') 'wprcmt 1 = ', wprcmt(1,1:5,is)
    write(*,'(1x,A,5ES18.10)') 'wprcmt 2 = ', wprcmt(2,1:5,is)
    write(*,'(1x,A,5ES18.10)') 'wprcmt 3 = ', wprcmt(3,1:5,is)
    write(*,'(1x,A,5ES18.10)') 'wprcmt 4 = ', wprcmt(4,1:5,is)

    write(*,*)
    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprcmt(2,1:5,is)
    !write(*,*)
    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprcmt(3,1:5,is)
    !write(*,*)
    !write(*,'(1x,A,5ES18.10)') 'wprmt = ', wprcmt(4,1:5,is)
  enddo

  call elk_stop()
end program
