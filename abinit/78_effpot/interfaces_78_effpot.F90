!!****m* ABINIT/interfaces_78_effpot
!! NAME
!! interfaces_78_effpot
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/78_effpot
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

module interfaces_78_effpot

 implicit none

interface
 subroutine compute_anharmonics(eff_pot,filenames,inp,comm)
  use m_multibinit_dataset
  use defs_basis
  use m_effective_potential
  implicit none
  integer, intent(in) :: comm
  type(effective_potential_type),target, intent(inout) :: eff_pot
  type(multibinit_dataset_type),intent(in) :: inp
  character(len=fnlen),intent(in) :: filenames(17)
 end subroutine compute_anharmonics
end interface

interface
 subroutine init10(filnam,comm)
  implicit none
  integer,intent(in) :: comm
  character(len=*),intent(out) :: filnam(17)
 end subroutine init10
end interface

end module interfaces_78_effpot
!!***
