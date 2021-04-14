! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlfv(nmatp,ngp,igpig,vgpc,apwalm,h)
  use modmain
  use modomp
  implicit none
  ! arguments
  integer, intent(in) :: nmatp,ngp,igpig(ngkmax)
  real(8), intent(in) :: vgpc(3,ngkmax)
  complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  complex(8), intent(out) :: h(nmatp,nmatp)
  ! local variables
  integer ias,i
  
  ! zero the upper triangular part of the matrix
  DO i=1,nmatp
    h(1:i,i) = 0.d0
  ENDDO 
  
  DO ias=1,natmtot
    CALL hmlaa(tefvr,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  ENDDO 
  CALL hmlistl(ngp,igpig,vgpc,nmatp,h)
  
  DO ias=1,natmtot
    CALL hmlalo(ias,ngp,apwalm(:,:,:,ias),nmatp,h)
    CALL hmllolo(ias,ngp,nmatp,h)
  ENDDO 
  
  RETURN 
END SUBROUTINE 

