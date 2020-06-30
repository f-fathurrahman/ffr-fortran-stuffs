include 'elk_init.f90'
include 'elk_stop.f90'
include 'prepare_all.f90'

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
