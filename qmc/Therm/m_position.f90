module m_position
! updates positions, their differences and two-particle potential and its change
  use m_random
  implicit none
  public :: JASEXP,DENSITY,DENSITY1
  integer,parameter,public :: NEMAX=100,NDIVMX=101
  integer,public :: NE,IE,NDIV,IMC
  real(kind=dp),parameter,public:: EMACH=1.0d-4
  real(kind=dp),public:: LENGTH
  real(kind=dp),dimension(NEMAX),public:: VNEW,VDIF
  real(kind=dp),dimension(3,NEMAX),public :: RE,RNEU
  real(kind=dp),dimension(4,NEMAX,NEMAX),public :: DIST,DISTNEU
  real(kind=dp),dimension(NDIVMX,NDIVMX,NDIVMX) :: RHO

contains

subroutine JASEXP(vn,vd)
!  updates the distances from the active electron to all others
        real(kind=dp),intent(out):: vn,vd
        integer :: i,k,n
        real(dp) :: work01,work02,work1,work2
        work1=0.0_dp
        work2=0.0_dp
        IENEK:do k=1,NE
         if (k .eq. IE) then
          cycle
         end if
         DIST(1:3,IE,k)=RE(1:3,IE)-RE(1:3,k)
         DISTNEU(1:3,IE,k)=RNEU(1:3,IE)-RE(1:3,k)
         work01=dsqrt(sum (DIST(1:3,IE,k)**2))
         work02=dsqrt(sum (DISTNEU(1:3,IE,k)**2))
         if (work01 .lt. EMACH)  work01 = 2.D0/3.D0*EMACH
         if (work02 .lt. EMACH)  work02 = 2.D0/3.D0*EMACH
         DIST(4,IE,k) = work01
         DISTNEU(4,IE,k) = work02
         work1=work1+1.D0/work02*dexp(-work02)-1.D0/work01*dexp(-work01)
         work2=work2+1.D0/work02*dexp(-work02)
        end do IENEK
        vd=work1
        vn=work2
end subroutine JASEXP

subroutine DENSITY(rh)
!  calculate the average particle density rh() on a cubic mesh
!  with NDIV intervals on each cubic axis, NDIV must be odd
        real(dp),intent(out),dimension(NDIVMX,NDIVMX,NDIVMX) :: rh
        integer :: nx,ny,nz,i
        real(dp) :: dl
        if (dble((NDIV-1)/2) .ne. dble(NDIV-1)/2.0_dp) then
         write(*,*) 'NDIV not odd: stop'
         stop
        end if
        rh = 0.0_dp
        dl = LENGTH/dble(NDIV-1)
        do i=1,NE
         nx = 1 + (NDIV-1)/2 + int(RE(1,i)/dl)
         ny = 1 + (NDIV-1)/2 + int(RE(2,i)/dl)
         nz = 1 + (NDIV-1)/2 + int(RE(3,i)/dl)
         if ((nx > NDIV) .or. (ny > NDIV) .or. (nz > NDIV)) then
          write(*,*)'too large nx, ny, or nz ',nx,ny,nz,' > NDIV ',NDIV
          stop
         end if
         rh(nx,ny,nz) = rh(nx,ny,nz) + 1
        end do
       end subroutine DENSITY

subroutine DENSITY1(NTOT,NDIV,LENGTH,RHO,NEMAX,NE,RE)
!  calculate the average particle density RHO() on a cubic mesh
!  with NDIV intervals on each cubic axis, NDIV must be even
       implicit none
       integer,intent(in)     :: NTOT,NDIV,NEMAX,NE
       real(dp),intent(in)    :: LENGTH,RE(3,NEMAX)
       real(dp),intent(inout) :: RHO(NTOT,NTOT,NTOT)
!  local variables
       integer  :: nx,ny,nz,i
       real(dp) :: dl
       if (DBLE(NDIV/2) .ne. DBLE(NDIV)/2.D0) THEN
        write(*,*) 'NDIV not even: stop'
        stop
       end if
       dl = LENGTH/DBLE(NDIV)
       do i=1,NE
         nx = 1 + NDIV/2 + INT(RE(1,i)/dl)
         ny = 1 + NDIV/2 + INT(RE(2,i)/dl)
         nz = 1 + NDIV/2 + INT(RE(3,i)/dl)
         if ((nx .gt. NTOT) .or. (ny .gt. NTOT) .or. (nz .gt. NTOT)) THEN
          write(*,*)'too large nx, ny, or nz ',nx,ny,nz,' > NTOT= ',NTOT
          stop
         end if
         RHO(nx,ny,nz) = RHO(nx,ny,nz) + 1
       end do
      end subroutine DENSITY1
end module
