!!****m* ABINIT/interfaces_28_numeric_noabirule
!! NAME
!! interfaces_28_numeric_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/28_numeric_noabirule
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_28_numeric_noabirule

 implicit none

interface
 subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )
  use defs_basis
  implicit none
  integer, intent(in) :: ntab
  real(dp), intent(in) :: a
  real(dp), intent(in) :: b
  real(dp), intent(out) :: result
  real(dp), intent(inout) :: e(ntab)
  real(dp), intent(in) :: ftab(ntab)
  real(dp), intent(inout) :: work(ntab)
  real(dp), intent(in) :: xtab(ntab)
  real(dp), intent(inout) :: y(3,ntab)
 end subroutine cspint
end interface

interface
 subroutine dzgedi(a,lda,n,ipvt,det,work,job)
  implicit none
  integer :: job
  integer :: lda
  integer :: n
  real*8 :: det(2,2)
  real*8 :: a(2,lda,n)
  integer :: ipvt(n)
  real*8 :: work(2,n)
 end subroutine dzgedi
end interface

interface
 subroutine dzgefa(a,lda,n,ipvt,info)
  implicit none
  integer :: info
  integer :: lda
  integer :: n
  real*8 :: a(2,lda,n)
  integer :: ipvt(n)
 end subroutine dzgefa
end interface

interface
 subroutine GAMMA_FUNCTION(X,GA)
  use defs_basis
  implicit none
  real(dp),intent(out) :: ga
  real(dp),intent(in) :: x
 end subroutine GAMMA_FUNCTION
end interface

interface
 function interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function interp
end interface

interface
 function dinterp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: dinterp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function dinterp
end interface

interface
 function taylor_interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: taylor_interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function taylor_interp
end interface

interface
 function dtaylor_interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: dtaylor_interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function dtaylor_interp
end interface

interface
 subroutine calculate_taylor_c(n,z,f,z0,c)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: z0
  complex(gwpc) :: c(n)
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end subroutine calculate_taylor_c
end interface

interface
 SUBROUTINE INTRPL(L,X,Y,N,U,V,dv,dv2,ideriv)
  implicit none
  integer, parameter :: NQQ=12000
  integer :: L
  integer :: N
  integer :: ideriv
  double precision :: DV(NQQ)
  double precision :: DV2(NQQ)
  double precision :: U(N)
  double precision :: V(N)
  double precision :: X(L)
  double precision :: Y(L)
 end subroutine INTRPL
end interface

interface
 SUBROUTINE CALJY0(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALJY0
end interface

interface
 DOUBLE PRECISION FUNCTION BESJ0(X)
 implicit none
 double precision :: X
end function BESJ0
end interface

interface
 DOUBLE PRECISION FUNCTION BESY0(X)
 implicit none
 double precision :: X
end function BESY0
end interface

interface
 SUBROUTINE CALJY1(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALJY1
end interface

interface
 DOUBLE PRECISION FUNCTION BESJ1(X)
 implicit none
 double precision :: X
end function BESJ1
end interface

interface
 DOUBLE PRECISION FUNCTION BESY1(X)
 implicit none
 double precision :: X
end function BESY1
end interface

interface
 subroutine jacobi(a,n,np,d,v,nrot)
  implicit none
  integer :: n
  integer :: np
  integer :: nrot
  real*8 :: a(np,np)
  real*8 :: d(np)
  real*8 :: v(np,np)
 end subroutine jacobi
end interface

interface
 SUBROUTINE CALCK0(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALCK0
end interface

interface
 DOUBLE PRECISION FUNCTION BESK0(X)
 implicit none
 double precision :: X
end function BESK0
end interface

interface
 DOUBLE PRECISION FUNCTION BESEK0(X)
 implicit none
 double precision :: X
end function BESEK0
end interface

interface
 SUBROUTINE CALCK1(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALCK1
end interface

interface
 DOUBLE PRECISION FUNCTION BESK1(X)
 implicit none
 double precision :: X
end function BESK1
end interface

interface
 DOUBLE PRECISION  FUNCTION BESEK1(X)
 implicit none
 double precision :: X
end function BESEK1
end interface

interface
 SUBROUTINE nfourier(rindata,coutdata,iflag,Iwmax,L,Beta)
  implicit none
  integer :: Iwmax
  integer :: L
  integer :: iflag
  double precision :: Beta
  complex*16 :: coutdata(Iwmax+1)
  double precision :: rindata(L)
 end subroutine nfourier
end interface

interface
 SUBROUTINE nfourier2(rindata,coutdata,iflag,om,L,Beta)
  implicit none
  integer :: L
  integer :: iflag
  double precision :: Beta
  complex*16 :: coutdata
  real*8 :: om
  double precision :: rindata(L)
 end subroutine nfourier2
end interface

interface
 SUBROUTINE invfourier(cindata,routdata,Iwmax,L,iflag,beta)
  implicit none
  integer, intent(in) :: Iwmax
  integer, intent(in) :: L
  integer, intent(in) :: iflag
  double precision, intent(in) :: beta
  complex*16, intent(in) :: cindata(1:Iwmax)
  complex*16, intent(inout) :: routdata(1:L)
 end subroutine invfourier
end interface

interface
 SUBROUTINE ludcmp(a,n,np,indx,id,info)
  implicit none
  integer :: id
  integer :: info
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  integer :: indx(n)
 end subroutine ludcmp
end interface

interface
 SUBROUTINE lubksb(a,n,np,indx,b)
  implicit none
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  real*8 :: b(n)
  integer :: indx(n)
 end subroutine lubksb
end interface


interface
 function uniformrandom(seed) 
  implicit none
  integer :: seed
  double precision :: uniformrandom
 end function uniformrandom
end interface

end module interfaces_28_numeric_noabirule
!!***
