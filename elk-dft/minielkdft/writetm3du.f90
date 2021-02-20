
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writetm3du
! !INTERFACE:
subroutine writetm3du(fnum)
! !USES:
use modmain
use moddftu
use modtest
! !DESCRIPTION:
!   Decompose the density matrix and the DFT+$U$ Hartree-Fock energy in 3-index
!   tensor moments components, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias,lm,i
integer l,p,k,r,rmin,t
real(8) ehf,vh,vx,w2,t1
real(8) tm3ptot,tm3p,tm3p0
! allocatable arrays
real(8), allocatable :: ex(:,:,:,:),tm32(:,:,:,:)
complex(8), allocatable :: tm3(:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(writetm3du): spin-polarisation must be enabled")')
  write(*,*)
  stop
end if
if (iscl.le.1) then
  write(fnum,*)
  write(fnum,'("Tensor moment decomposition of density matrix and Hartree-Fock &
   &energy")')
  write(fnum,'(" (see Phys. Rev. B. 80, 035121 (2009))")')
  write(fnum,'("  W.W = modulus square of tensor moment")')
  write(fnum,'("  Eh = Hartree energy term")')
  write(fnum,'("  Ex = exchange energy term")')
  write(fnum,'("  Pol = polarisation of density matrix")')
end if
if (iscl.ge.1) then
  write(fnum,*)
  write(fnum,'("+--------------------+")')
  write(fnum,'("| Loop number : ",I4," |")') iscl
  write(fnum,'("+--------------------+")')
end if
! allocate arrays
lm=2*lmaxdm+1
allocate(ex(0:2*lmaxdm,0:1,0:lm,natmtot))
allocate(tm32(0:2*lmaxdm,0:1,0:lm,natmtot))
allocate(tm3(-lmmaxdm:lmmaxdm))
tm32(:,:,:,:)=0.d0
ex(:,:,:,:)=0.d0
tm3p0=0.d0
! loop over DFT+U entries
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! write info to TMDFTU.OUT
    write(fnum,*)
    write(fnum,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
     trim(spsymb(is)),ia
    write(fnum,'(" l = ",I2)') l
! zero total Hartree-Fock energy
    ehf=0.d0
! zero total polarisation
    tm3ptot=0.d0
    write(fnum,*)
    do k=0,2*l
      do p=0,1
        rmin=abs(k-p)
        do r=rmin,k+p
! decompose density matrix in 3-index tensor moment components
          call dmtotm3(l,nspinor,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),tm3)
          w2=0.d0
! calculate modulus square of tensor moment
          do t=-r,r
            w2=w2+dble((-1)**t)*dble(tm3(-t)*tm3(t))
          end do
! potential energy components
          call pottm3(i,k,p,r,vh,vx)
! save square of tensor modulus and exchange energy
          tm32(k,p,r,ias)=w2
          ex(k,p,r,ias)=vx*w2
! polarisation terms
          call tm3pol(l,k,p,r,w2,tm3p)
! write square of tensor modulus; Hartree and exchange energy; and polarisation
          write(fnum,'("  k = ",I2,", p = ",I2,", r = ",I2)') k,p,r
          if ((k+p+r).eq.0) then
! for k,p,r = 0 save reference polarisation
            tm3p0=tm3p
! for k,p,r = 0 do not write out the polarisation
            write(fnum,'("   W.W =",F14.8,", Eh =",F14.8,", Ex =",F14.8)') w2, &
             vh*w2,vx*w2
          else
! relative polarisation
            t1=tm3p/tm3p0
            write(fnum,'("   W.W =",F14.8,", Eh =",F14.8,", Ex =",F14.8,&
             &", Pol =",F14.8)') w2,vh*w2,vx*w2,t1
! total relative polarisation (skipping 000 component)
            tm3ptot=tm3ptot+t1
          end if
          ehf=ehf+(vh+vx)*w2
          do t=-r,r
! write out single components of tensor moments
            write(fnum,'("    t = ",I2," : ",2F14.8)') t,tm3(t)
          end do
          write(fnum,*)
        end do
      end do
    end do
    write(fnum,*)
    write(fnum,'("  Total Hartree-Fock energy (without DC correction) : ",&
     &F14.8)') ehf
    write(fnum,'("  Total polarisation of density matrix : ",F14.8)') tm3ptot
    write(fnum,*)
! end loop over atoms and species
  end do
end do
flush(fnum)
! write test files if required
if (test) then
  t1=sqrt(sum(tm32(:,:,:,:)**2))
  call writetest(820,'Norm of tensor moments',tol=1.d-4,rv=t1)
  t1=sqrt(sum(ex(:,:,:,:)**2))
  call writetest(830,'Norm of DFT+U Hartree-Fock exchange energies',tol=1.d-4, &
   rv=t1)
end if
deallocate(ex,tm32,tm3)
return
end subroutine
!EOC

