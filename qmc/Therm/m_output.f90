module m_output
!  output for 1-,2-,3-dimensional arrays on files termed PRONAME
!  with some additional naming
  use m_random
  implicit none
  public :: OUT3D
  integer,parameter,public :: FE1MAX=100,FE2MAX=100,FE3MAX=100
  integer :: ios
  character(7),public :: PRONAME

contains

subroutine OUT3D(NDIV1,NDIV2,NDIV3,a3d)
  integer :: NDIV1,NDIV2,NDIV3
  real(kind=dp),intent(in) :: a3d(FE1MAX,FE2MAX,FE3MAX)
  integer :: nx,ny,nz
  open(unit=36,file=PRONAME//"_erw10.dat",position="append",status="unknown")
        do nz=1,NDIV3
         do ny=1,NDIV2
          do nx=1,NDIV1
           write(36,55) nx,ny,nz,a3d(nx,ny,nz)
          end do
         end do
         write(36,*)'  '
         write(36,*)'  '
        end do
        close(36,iostat=ios)
        if (ios .gt. 0) then
          write(35,*) 'output error from OUT3D'
        end if
55     format(t3,i3,i3,i3,e12.3)
       end subroutine OUT3D
end module
