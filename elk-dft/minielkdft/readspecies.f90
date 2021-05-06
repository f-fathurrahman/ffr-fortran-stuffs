
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

SUBROUTINE readspecies

!use modmain
use modmain, only: sppath, spfname, spsymb, spzn, spmass, nspecies, &
                   spname, spcore, rminsp, nstsp, nsp, ksp, lsp, occsp, &
                   nrmt, rmt, rmaxsp, nrmtscf, maxstsp, nxoapwlo, nxlo, rmtall, &
                   apword, apwe0, apwdm, apwve, nlorb, lorbord, lorbve, lorbdm, &
                   lorbe0, lorbl, maxlorbord, maxlorb, lmaxo, maxapword, lmaxapw, e0min
IMPLICIT NONE 
! local variables
INTEGER is,ist,ios
INTEGER nlx,ilx,lx,ilo
INTEGER io,jo,ko,l,i,j

e0min=0.d0
DO is=1,nspecies
  open(50,file=trim(sppath)//trim(spfname(is)),status='OLD',form='FORMATTED', &
   iostat=ios)
  IF(ios.ne.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): error opening species file ",A)') &
     trim(sppath)//trim(spfname(is))
    WRITE(*,*)
    stop
  ENDIF 
  read(50,*) spsymb(is)
  read(50,*) spname(is)
  read(50,*) spzn(is)
  read(50,*) spmass(is)
  read(50,*) rminsp(is),rmt(is),rmaxsp(is),nrmt(is)
  IF(rminsp(is).le.0.d0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): rminsp <= 0 : ",G18.10)') rminsp(is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 
  IF(rmt(is).le.rminsp(is)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): rmt <= rminsp : ",2G18.10)') rmt(is), &
     rminsp(is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 
  IF(rmaxsp(is).lt.rmt(is)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): rmaxsp < rmt : ",2G18.10)') rmaxsp(is), &
     rmt(is)
    WRITE(*,*)
    stop
  ENDIF 
  IF(nrmt(is).lt.20) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): nrmt too small : ",I8)') nrmt(is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 
! multiply nrmt by the scale factor
  nrmt(is)=nint(dble(nrmt(is))*nrmtscf)
! reduce the minimum radial mesh point by the same factor
  rminsp(is)=rminsp(is)/nrmtscf
  read(50,*) nstsp(is)
  IF((nstsp(is).le.0).or.(nstsp(is).gt.maxstsp)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): nstsp out of range : ",I8)') nstsp(is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 
  DO ist=1,nstsp(is)
    read(50,*) nsp(ist,is),lsp(ist,is),ksp(ist,is),occsp(ist,is),spcore(ist,is)
    IF(nsp(ist,is).lt.1) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): nsp < 1 : ",I8)') nsp(ist,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and state ",I4)') ist
      WRITE(*,*)
      stop
    ENDIF 
    IF(lsp(ist,is).lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lsp < 0 : ",I8)') lsp(ist,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and state ",I4)') ist
      WRITE(*,*)
      stop
    ENDIF 
    IF(ksp(ist,is).lt.1) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): ksp < 1 : ",I8)') ksp(ist,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and state ",I4)') ist
      WRITE(*,*)
      stop
    ENDIF 
    IF(occsp(ist,is).lt.0.d0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): occsp < 0 : ",G18.10)') occsp(ist,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and state ",I4)') ist
      WRITE(*,*)
      stop
    ENDIF 
  ENDDO 

  read(50,*) apword(0,is)
  IF(apword(0,is).le.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(0,is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 
  IF(apword(0,is).gt.maxapword) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): apword too large : ",I8)') apword(0,is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,'("Adjust maxapword in modmain and recompile code")')
    WRITE(*,*)
    stop
  ENDIF 
! set the APW orders for l>0
  apword(1:lmaxapw,is)=apword(0,is)
  DO io=1,apword(0,is)
    read(50,*) apwe0(io,0,is),apwdm(io,0,is),apwve(io,0,is)
    IF(apwdm(io,0,is).lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,0,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and order ",I4)') io
      WRITE(*,*)
      stop
    ENDIF 
! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io,1:lmaxapw,is)=apwe0(io,0,is)
    apwdm(io,1:lmaxapw,is)=apwdm(io,0,is)
    apwve(io,1:lmaxapw,is)=apwve(io,0,is)
    e0min=min(e0min,apwe0(io,0,is))
  ENDDO 

  read(50,*) nlx
  IF(nlx.lt.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): nlx < 0 : ",I8)') nlx
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 

  ! Read exceptions if any
  DO ilx=1,nlx
    read(50,*) lx,io
    IF(lx.lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lx < 0 : ",I8)') lx
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and exception number ",I4)') ilx
      WRITE(*,*)
      stop
    ENDIF 
    IF(lx.gt.lmaxapw) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lx > lmaxapw : ",I8)') lx
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and exception number ",I4)') ilx
      WRITE(*,*)
      stop
    ENDIF 
    apword(lx,is)=io
    IF(apword(lx,is).le.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(lx,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and exception number ",I4)') ilx
      WRITE(*,*)
      stop
    ENDIF 
    IF(apword(lx,is).gt.maxapword) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): apword too large : ",I8)') apword(lx,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and exception number ",I4)') ilx
      WRITE(*,'("Adjust maxapword in modmain and recompile code")')
      WRITE(*,*)
      stop
    ENDIF 
    DO io=1,apword(lx,is)
      read(50,*) apwe0(io,lx,is),apwdm(io,lx,is),apwve(io,lx,is)
      IF(apwdm(io,lx,is).lt.0) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,lx,is)
        WRITE(*,'(" for species ",I4)') is
        WRITE(*,'(" exception number ",I4)') ilx
        WRITE(*,'(" and order ",I4)') io
        WRITE(*,*)
        stop
      ENDIF 
      e0min=min(e0min,apwe0(io,lx,is))
    ENDDO 
  ENDDO 

! add excess order to APW functions if required
  IF(nxoapwlo.gt.0) THEN 
    DO l=0,lmaxapw
      jo=apword(l,is)
      ko=jo+nxoapwlo
      IF(ko.gt.maxapword) ko=maxapword
      i=0
      DO io=jo+1,ko
        i=i+1
        apwe0(io,l,is)=apwe0(jo,l,is)
        apwdm(io,l,is)=apwdm(jo,l,is)+i
        apwve(io,l,is)=apwve(jo,l,is)
      ENDDO 
      apword(l,is)=ko
    ENDDO 
  ENDIF 

  read(50,*) nlorb(is)
  IF(nlorb(is).lt.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): nlorb < 0 : ",I8)') nlorb(is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,*)
    stop
  ENDIF 
  IF(nlorb(is).gt.maxlorb) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(readspecies): nlorb too large : ",I8)') nlorb(is)
    WRITE(*,'(" for species ",I4)') is
    WRITE(*,'("Adjust maxlorb in modmain and recompile code")')
    WRITE(*,*)
    stop
  ENDIF 
  DO ilo=1,nlorb(is)
    read(50,*) lorbl(ilo,is),lorbord(ilo,is)
    IF(lorbl(ilo,is).lt.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lorbl < 0 : ",I8)') lorbl(ilo,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and local-orbital ",I4)') ilo
      WRITE(*,*)
      stop
    ENDIF 
    IF(lorbl(ilo,is).gt.lmaxo) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lorbl > lmaxo : ",2I8)') lorbl(ilo,is), &
       lmaxo
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and local-orbital ",I4)') ilo
      WRITE(*,*)
      stop
    ENDIF 
    IF(lorbord(ilo,is).lt.2) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lorbord < 2 : ",I8)') lorbord(ilo,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and local-orbital ",I4)') ilo
      WRITE(*,*)
      stop
    ENDIF 
    IF(lorbord(ilo,is).gt.maxlorbord) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(readspecies): lorbord too large : ",I8)') lorbord(ilo,is)
      WRITE(*,'(" for species ",I4)') is
      WRITE(*,'(" and local-orbital ",I4)') ilo
      WRITE(*,'("Adjust maxlorbord in modmain and recompile code")')
      WRITE(*,*)
      stop
    ENDIF 
    DO io=1,lorbord(ilo,is)
      read(50,*) lorbe0(io,ilo,is),lorbdm(io,ilo,is),lorbve(io,ilo,is)
      IF(lorbdm(io,ilo,is).lt.0) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(readspecies): lorbdm < 0 : ",I8)') lorbdm(io,ilo,is)
        WRITE(*,'(" for species ",I4)') is
        WRITE(*,'(" local-orbital ",I4)') ilo
        WRITE(*,'(" and order ",I4)') io
        WRITE(*,*)
        stop
      ENDIF 
      e0min=min(e0min,lorbe0(io,ilo,is))
    ENDDO 
  ENDDO 
! add excess local-orbitals if required
  IF(nxlo.gt.0) THEN 
    lx=-1
    DO ilo=1,nlorb(is)
      DO io=1,lorbord(ilo,is)
        IF(lorbe0(io,ilo,is).lt.0.d0) goto 10
      ENDDO 
      IF(lorbl(ilo,is).gt.lx) lx=lorbl(ilo,is)
10 continue
    ENDDO 
    ilo=nlorb(is)
    DO i=1,nxlo
      IF(ilo.eq.maxlorb) exit
      l=lx+i
      IF(l.gt.lmaxo) exit
      ilo=ilo+1
      lorbl(ilo,is)=l
      lorbord(ilo,is)=apword(l,is)+1
      DO io=1,lorbord(ilo,is)
        lorbe0(io,ilo,is)=apwe0(1,l,is)
        lorbdm(io,ilo,is)=io-1
        lorbve(io,ilo,is)=apwve(1,l,is)
      ENDDO 
    ENDDO 
    nlorb(is)=ilo
  ENDIF 
! add excess order to local-orbitals if required
  IF(nxoapwlo.gt.0) THEN 
    DO ilo=1,nlorb(is)
! find the maximum energy derivative
      jo=1
      j=lorbdm(jo,ilo,is)
      DO io=1,lorbord(ilo,is)
        i=lorbdm(io,ilo,is)
        IF(i.gt.j) THEN 
          jo=io
          j=i
        ENDIF 
      ENDDO 
      ko=lorbord(ilo,is)+nxoapwlo
      IF(ko.gt.maxlorbord) ko=maxlorbord
      i=0
      DO io=lorbord(ilo,is)+1,ko
        i=i+1
        lorbe0(io,ilo,is)=lorbe0(jo,ilo,is)
        lorbdm(io,ilo,is)=lorbdm(jo,ilo,is)+i
        lorbve(io,ilo,is)=lorbve(jo,ilo,is)
      ENDDO 
      lorbord(ilo,is)=ko
    ENDDO 
  ENDIF 
  close(50)
ENDDO 

! set all muffin-tin radii to single value if required
IF(rmtall.gt.0.d0) rmt(1:nspecies)=rmtall

! add conduction state local-orbitals if required
CALL addlorbcnd

! subtract 2 Hartree from the minimum energy
e0min=e0min-2.d0

RETURN 
END SUBROUTINE 

