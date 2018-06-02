!!****m* ABINIT/interfaces_61_occeig
!! NAME
!! interfaces_61_occeig
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/61_occeig
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

module interfaces_61_occeig

 implicit none

interface
 subroutine dos_hdr_write(buffer,deltaene,dosdeltae,&  
  &  eigen,enemax,enemin,fermie,mband,nband,nene,&  
  &  nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: prtdos
  integer,intent(in) :: unitdos
  real(dp),intent(in) :: buffer
  real(dp),intent(in) :: deltaene
  real(dp),intent(in) :: dosdeltae
  real(dp),intent(in) :: enemax
  real(dp),intent(in) :: enemin
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine dos_hdr_write
end interface

interface
 subroutine getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&  
  &  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: unitdos
  real(dp),intent(in) :: dosdeltae
  real(dp),intent(out) :: entropy
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: maxocc
  real(dp),intent(out) :: nelect
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine getnel
end interface

interface
 subroutine init_occ_ent(entfun,limit,nptsdiv2,occfun,occopt,option,smdfun,tphysel,tsmear,tsmearinv,xgrid)
  use defs_basis
  implicit none
  integer,intent(inout) :: nptsdiv2
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  real(dp),intent(out) :: limit
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: tsmearinv
  real(dp),intent(inout) :: entfun(-nptsdiv2:nptsdiv2,2)
  real(dp),intent(inout) :: occfun(-nptsdiv2:nptsdiv2,2)
  real(dp),intent(inout) :: smdfun(-nptsdiv2:nptsdiv2,2)
  real(dp),intent(inout) :: xgrid(-nptsdiv2:nptsdiv2)
 end subroutine init_occ_ent
end interface

interface
 subroutine newocc(doccde,eigen,entropy,fermie,spinmagntarget,mband,nband,&  
  &  nelect,nkpt,nspinor,nsppol,occ,occopt,prtvol,stmbias,tphysel,tsmear,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  real(dp),intent(in) :: nelect
  real(dp),intent(in) :: spinmagntarget
  real(dp),intent(in) :: stmbias
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine newocc
end interface

interface
 subroutine occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,occopt,occ_k,occ_kq,rocceig)
  use defs_basis
  implicit none
  integer,intent(in) :: nband_k
  integer,intent(in) :: occopt
  real(dp),intent(in) :: doccde_k(nband_k)
  real(dp),intent(in) :: doccde_kq(nband_k)
  real(dp),intent(in) :: eig0_k(nband_k)
  real(dp),intent(in) :: eig0_kq(nband_k)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: occ_kq(nband_k)
  real(dp),intent(out) :: rocceig(nband_k,nband_k)
 end subroutine occeig
end interface

interface
 subroutine pareigocc(eigen,formeig,localrdwf,mpi_enreg,mband,nband,nkpt,nsppol,occ,transmit_occ)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: localrdwf
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: transmit_occ
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(inout) :: eigen(mband*(2*mband)**formeig*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 end subroutine pareigocc
end interface

end module interfaces_61_occeig
!!***
