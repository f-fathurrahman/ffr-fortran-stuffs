!-------------------------------
SUBROUTINE init_gvector_arrays()
!-------------------------------

  USE modmain, ONLY: &
               nspecies, isgkmax, rmt, vgc, sfacg, rgkmax, npsd, natoms, ngvec, &
               natmtot, lnpsd, lmaxo, gmaxvr, gkmax, jlgrmt, gc, ffacg, ngtot

  IMPLICIT NONE 
  INTEGER :: is
  REAL(8) :: t1, rsum

  !-------------------------!
  !     G-vector arrays     !
  !-------------------------!
  ! determine gkmax from rgkmax
  IF (nspecies.eq.0) isgkmax=-2
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
  CALL gengvec()
  
  ! write number of G-vectors to test file
  WRITE(*,*) 'number of G-vectors = ', ngvec
  
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
  !call writevars('gmaxvr',rv=gmaxvr)
  !call writevars('ngridg',nv=3,iva=ngridg)
  !call writevars('intgv',nv=6,iva=intgv)
  !call writevars('ngvec',iv=ngvec)
  !call writevars('ivg',nv=3*ngtot,iva=ivg)
  !call writevars('igfft',nv=ngtot,iva=igfft)

END SUBROUTINE 

