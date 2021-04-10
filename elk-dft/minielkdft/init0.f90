SUBROUTINE init0()
  USE modmain
  USE modxcifc
  USE moddftu
  USE modtddft
  USE modphonon
  USE modulr
  USE modtest
  USE modvars
  USE modmpi
  USE modomp
  !   Performs basic consistency checks as well as allocating and initialising
  !   global variables not dependent on the $k$-point set.
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ias,ist
  INTEGER :: nr,l,i
  REAL(8) :: rsum,t1
  REAL(8) :: ts0,ts1
  
  WRITE(*,*)
  WRITE(*,*) '*** ffr: Entering init0 ***'
  WRITE(*,*)
  
  
  !-------------------------------!
  !     zero timing variables     !
  !-------------------------------!
  timeinit=0.d0
  timemat=0.d0
  timefv=0.d0
  timesv=0.d0
  timerho=0.d0
  timepot=0.d0
  timefor=0.d0
  CALL timesec(ts0)
  
  CALL init_am_variables()
  
  CALL init_idx_atom_species()
  
  CALL init_spin_variables()

  CALL init_crystal_structure()

  CALL init_vector_field_E_A()

  !---------------------------------!
  !     crystal symmetry set up     !
  !---------------------------------!
  CALL symmetry()

  CALL init_radial_meshes()

!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
chgzn=0.d0
chgcrtot=0.d0
chgval=0.d0
nstspmax=0
nstcr=0
do is=1,nspecies
! nuclear charge
  chgzn=chgzn+spzn(is)*natoms(is)
! find the maximum number of atomic states
  nstspmax=max(nstspmax,nstsp(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
  spze(is)=0.d0
  chgcr(is)=0.d0
  do ist=1,nstsp(is)
    spze(is)=spze(is)+occsp(ist,is)
    if (spcore(ist,is)) then
      chgcr(is)=chgcr(is)+occsp(ist,is)
      nstcr=nstcr+2*ksp(ist,is)*natoms(is)
    else
      chgval=chgval+occsp(ist,is)*natoms(is)
    end if
  end do
  chgcrtot=chgcrtot+chgcr(is)*natoms(is)
end do
! add excess charge
chgval=chgval+chgexs
! total charge
chgtot=chgcrtot+chgval
if (chgtot.lt.1.d-8) then
  write(*,*)
  write(*,'("Error(init0): zero total charge")')
  write(*,*)
  stop
end if
! effective Wigner radius
rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)
! write to VARIABLES.OUT
call writevars('spze',nv=nspecies,rva=spze)
call writevars('chgcr',nv=nspecies,rva=chgcr)
call writevars('chgexs',rv=chgexs)
call writevars('chgval',rv=chgtot)

!-------------------------!
!     G-vector arrays     !
!-------------------------!
! determine gkmax from rgkmax
if (nspecies.eq.0) isgkmax=-2
select case(isgkmax)
case(:-4)
! use largest muffin-tin radius
  gkmax=rgkmax/maxval(rmt(1:nspecies))
case(-3)
! use smallest muffin-tin radius
  gkmax=rgkmax/minval(rmt(1:nspecies))
case(-2)
! use the fixed value of 2.0
  gkmax=rgkmax/2.d0
case(-1)
! use average muffin-tin radius
  rsum=0.d0
  do is=1,nspecies
     rsum=rsum+dble(natoms(is))*rmt(is)
  end do
  rsum=rsum/dble(natmtot)
  gkmax=rgkmax/rsum
case(1:)
! use user-specified muffin-tin radius
  if (isgkmax.le.nspecies) then
    gkmax=rgkmax/rmt(isgkmax)
  else
    write(*,*)
    write(*,'("Error(init0): isgkmax > nspecies : ",2I8)') isgkmax,nspecies
    write(*,*)
    stop
  end if
end select
! generate the G-vectors
call gengvec
! write number of G-vectors to test file
call writetest(900,'number of G-vectors',iv=ngvec)
! Poisson solver pseudocharge density constant
if (nspecies.gt.0) then
  t1=0.25d0*gmaxvr*maxval(rmt(1:nspecies))
else
  t1=0.25d0*gmaxvr*2.d0
end if
npsd=max(nint(t1),1)
lnpsd=lmaxo+npsd+1

! generate the Coulomb Green's function in G-space = fourpi / G^2
call gengclg()

! compute the spherical Bessel functions j_l(|G|R_mt)
if (allocated(jlgrmt)) deallocate(jlgrmt)
allocate(jlgrmt(0:lnpsd,ngvec,nspecies))
call genjlgprmt(lnpsd,ngvec,gc,ngvec,jlgrmt)

! generate the spherical harmonics of the G-vectors
call genylmg()

! allocate structure factor array for G-vectors
if (allocated(sfacg)) deallocate(sfacg)
allocate(sfacg(ngvec,natmtot))

! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)

! generate the smooth step function form factors
if (allocated(ffacg)) deallocate(ffacg)
allocate(ffacg(ngtot,nspecies))
do is=1,nspecies
  call genffacgp(is,gc,ffacg(:,is))
end do

! generate the characteristic function
call gencfun()

! write to VARIABLES.OUT
call writevars('gmaxvr',rv=gmaxvr)
call writevars('ngridg',nv=3,iva=ngridg)
call writevars('intgv',nv=6,iva=intgv)
call writevars('ngvec',iv=ngvec)
call writevars('ivg',nv=3*ngtot,iva=ivg)
call writevars('igfft',nv=ngtot,iva=igfft)

!-------------------------!
!     atoms and cores     !
!-------------------------!

! determine the nuclear Coulomb potential
if (allocated(vcln)) deallocate(vcln)
allocate(vcln(nrspmax,nspecies))
t1=1.d0/y00
do is=1,nspecies
  nr=nrsp(is)
  call potnucl(ptnucl,nr,rsp(:,is),spzn(is),vcln(:,is))
  vcln(1:nr,is)=t1*vcln(1:nr,is)
end do

! solve the Kohn-Sham-Dirac equations for all atoms
call allatoms()

! allocate core state occupancy and eigenvalue arrays and set to default
if (allocated(occcr)) deallocate(occcr)
allocate(occcr(nstspmax,natmtot))
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(nstspmax,natmtot))
do ias=1,natmtot
  is=idxis(ias)
  do ist=1,nstsp(is)
    occcr(ist,ias)=occsp(ist,is)
    evalcr(ist,ias)=evalsp(ist,is)
  end do
end do

! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(nrspmax,2,nstspmax,natmtot))

! number of core spin channels
if (spincore) then
  nspncr=2
else
  nspncr=1
end if

! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(nrmtmax,natmtot,nspncr))

!-------------------------------------------------------------!
!     charge density, potentials and exchange-correlation     !
!-------------------------------------------------------------!
! allocate charge density arrays
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(npmtmax,natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngtot))
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (spinpol) then
  allocate(magmt(npmtmax,natmtot,ndmag))
  allocate(magir(ngtot,ndmag))
end if
! check if the current density should be calculated
if (tafield.or.(any(task.eq.[371,372,373,460,461]))) then
  tcden=.true.
end if
! allocate current density arrays
if (allocated(cdmt)) deallocate(cdmt)
if (allocated(cdir)) deallocate(cdir)
if (tcden) then
  allocate(cdmt(npmtmax,natmtot,3),cdir(ngtot,3))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(npmtmax,natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngtot))
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(npmtmax,natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(npmtmax,natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(npmtmax,natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngtot))
! effective Kohn-Sham potential
if (allocated(vsmt)) deallocate(vsmt)
allocate(vsmt(npmtmax,natmtot))
if (allocated(vsir)) deallocate(vsir)
allocate(vsir(ngtot))
if (allocated(vsig)) deallocate(vsig)
allocate(vsig(ngvec))
! exchange-correlation, dipole and Kohn-Sham effective magnetic fields
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (allocated(bdmt)) deallocate(bdmt)
if (allocated(bdir)) deallocate(bdir)
if (allocated(bsmt)) deallocate(bsmt)
if (allocated(bsir)) deallocate(bsir)
if (spinpol) then
  allocate(bxcmt(npmtmax,natmtot,ndmag),bxcir(ngtot,ndmag))
  if (tbdip) allocate(bdmt(npmtmax,natmtot,ndmag),bdir(ngtot,ndmag))
  allocate(bsmt(npcmtmax,natmtot,ndmag),bsir(ngtot,ndmag))
end if
! kinetic energy density
if (allocated(taumt)) deallocate(taumt)
if (allocated(tauir)) deallocate(tauir)
if (allocated(taucr)) deallocate(taucr)
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) then
  allocate(taumt(npmtmax,natmtot,nspinor),tauir(ngtot,nspinor))
  allocate(taucr(npmtmax,natmtot,nspinor))
end if
! tau-DFT exchange-correlation and Kohn-Sham potentials
if (allocated(wxcmt)) deallocate(wxcmt)
if (allocated(wxcir)) deallocate(wxcir)
if (allocated(wsmt)) deallocate(wsmt)
if (allocated(wsir)) deallocate(wsir)
if (xcgrad.eq.4) then
  allocate(wxcmt(npmtmax,natmtot),wxcir(ngtot))
  allocate(wsmt(npcmtmax,natmtot),wsir(ngtot))
  tevecsv=.true.
end if
! spin-orbit coupling radial function
if (allocated(socfr)) deallocate(socfr)
if (spinorb) then
  allocate(socfr(nrcmtmax,natmtot))
end if
! allocate muffin-tin charge and moment arrays
if (allocated(chgcrlk)) deallocate(chgcrlk)
allocate(chgcrlk(natmtot))
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3,natmtot))
! check if scaled spin exchange-correlation (SSXC) should be used
if (abs(ssxc-1.d0).gt.1.d-6) then
  tssxc=.true.
else
  tssxc=.false.
end if

!-------------------------!
!     force variables     !
!-------------------------!
if (allocated(forcehf)) deallocate(forcehf)
allocate(forcehf(3,natmtot))
if (allocated(forceibs)) deallocate(forceibs)
allocate(forceibs(3,natmtot))
if (allocated(forcetot)) deallocate(forcetot)
allocate(forcetot(3,natmtot))

!-------------------------------------------------!
!     DFT+U and fixed tensor moment variables     !
!-------------------------------------------------!
if ((dftu.ne.0).or.(ftmtype.ne.0)) then
! density matrix elements in each muffin-tin
  if (allocated(dmatmt)) deallocate(dmatmt)
  allocate(dmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! potential matrix elements in each muffin-tin
  if (allocated(vmatmt)) deallocate(vmatmt)
  allocate(vmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! zero the potential matrix
  vmatmt(:,:,:,:,:)=0.d0
! require the potential matrix elements be calculated
  tvmatmt=.true.
! flags for non-zero muffin-tin potential matrices
  if (allocated(tvmmt)) deallocate(tvmmt)
  allocate(tvmmt(0:lmaxdm,natmtot))
  tvmmt(:,:)=.false.
! require second-variational eigenvectors
  tevecsv=.true.
end if
if (dftu.ne.0) then
! DFT+U energy for each atom
  if (allocated(engyadu)) deallocate(engyadu)
  allocate(engyadu(natmmax,ndftu))
! interpolation constants (alpha)
  if (allocated(alphadu)) deallocate(alphadu)
  allocate(alphadu(natmmax,ndftu))
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ndftu
    is=idftu(1,i)
    if (is.gt.nspecies) then
      write(*,*)
      write(*,'("Error(init0): invalid species number : ",I8)') is
      write(*,*)
      stop
    end if
    l=idftu(2,i)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      tvmmt(l,ias)=.true.
    end do
  end do
end if
if (ftmtype.ne.0) then
! allocate and zero the fixed tensor moment potential array
  if (allocated(vmftm)) deallocate(vmftm)
  allocate(vmftm(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
  vmftm(:,:,:,:,:)=0.d0
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ntmfix
    is=itmfix(1,i)
    ia=itmfix(2,i)
    ias=idxas(ia,is)
    l=itmfix(3,i)
    tvmmt(l,ias)=.true.
  end do
end if

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine nuclear radii and volumes
call nuclei
! determine the nuclear-nuclear energy
call energynn
! get smearing function description
call getsdata(stype,sdescr)
! get mixing type description
call getmixdata(mixtype,mixdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat
! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! initial self-consistent loop number
iscl=1
tlast=.false.
! set the Fermi energy to zero
efermi=0.d0
! set the temperature from the smearing width
tempk=swidth/kboltz

call timesec(ts1)
timeinit=timeinit+ts1-ts0

write(*,*)
write(*,*) '*** ffr: Leaving init0 ***'
write(*,*)

return
end subroutine
!EOC

