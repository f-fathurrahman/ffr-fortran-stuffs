!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  A MC run for the Hydrogen atom
!  a.u. used
!
!-----------------------------------------------------------------------
!  A full QMC run for the H atom
program HMC
  use random
  use position
  implicit none
  integer, parameter    :: NXFA1=59,NXFA2=59
  real(dp)             :: STEPMAX,ALPHA,DALPHA,ALPHA0,ALPHA1
  real(dp)             :: WAVEC,DWAVEC,WAVEC0,WAVEC1
  real(dp)             :: RDIF,QUOT,LOCEN
  real(dp)             :: ERWEN,VAREN,ERWKIN,VARKIN,ERWPOT,VARPOT
  real(dp)             :: RAD1,RAD2,RADNEU,EK,EP
  logical              :: MCSCHRITT
  !  Local variables
  integer                     :: ix,n1,n2,n3
  real(dp)                    :: rannumb
  real(dp),dimension(NXFA2+1) :: work1,work2
  
  !   Controll output
  open(unit=35,FILE="MC21.OUT",STATUS="UNKNOWN")
  open(unit=38,FILE="Hmcvar_splot21.dat",STATUS="UNKNOWN")
  open(unit=39,FILE="Hmcerw_splot21.dat",STATUS="UNKNOWN")
  
  write(36,*) 'Variance on (c,alpha) - plane'
  write(36,*) NXFA2+1,'  ',NXFA1+1
  write(37,*) 'Energy expectation value on (c,alpha) - plane'
  write(37,*) NXFA2+1,'  ',NXFA1+1
  
  ! Number MC steps
  MCPRE = 1000000
  MCMAX = 1000000
  DRHO = 0.01_dp
  
  ! Start data:
  lnxfa1:do n1 = 1,NXFA1+1
  lnxfa2:do n2 = 1,NXFA2+1
    COUNTU = 0
    COUNTD = 0
    
    ! Wave function coefficients, 0<ALPHA<=2., -1.<=WAVEC<=2.
    ALPHA0=0.1_dp
    ALPHA1=1.9_dp
    WAVEC0=-0.8_dp
    WAVEC1=+1.0_dp
    DALPHA=(ALPHA1-ALPHA0)/DBLE(NXFA1+1)
    DWAVEC=(WAVEC1-WAVEC0)/DBLE(NXFA2+1)
    ALPHA=ALPHA0+(n1-1)*DALPHA
    WAVEC=WAVEC0 + (n2-1)*DWAVEC
    
    call INITRAN()
  
    ! Initial electron position
    RE(1) = 0.1_dp
    RE(2) = 0.1_dp
    RE(3) = 0.1_dp
    RNEU = RE
    
    ! Maximum step width, KONTUZ: check with acceptance ratio!
    STEPMAX = 2.D-0/(1.D0+ALPHA)
    
    ! Counts the acceptance number
    MCOUNT = 0
    
    ! Observables
    LOCEN = 0.0_dp
    ERWEN = 0.0_dp
    VAREN = 0.0_dp
    ERWKIN = 0.0_dp
    VARKIN = 0.0_dp
    ERWPOT = 0.0_dp
    VARPOT = 0.0_dp
    AVRHORAD = 0.0_dp


    ! MC loop: prerun for thermalizing,
    ! KONTUZ: does not change the bad sampling of energy!!!
    lrunpre:do IMC=1,MCPRE
      
      do ix=1,3
        ! Shift position at random within +-STEPMAX/2
        call GENRAN(rannumb)
        RDIF = (rannumb-0.5)*STEPMAX
        RNEU(ix) = RE(ix)+RDIF
      end do
      
      ! Calculate wave function ratio psi=(1+c*r)exp(-alpha*r)
      RAD1 = DSQRT(RE(1)**2+RE(2)**2+RE(3)**2)
      RAD2 = DSQRT(RNEU(1)**2+RNEU(2)**2+RNEU(3)**2)
      QUOT=((1.0_dp+WAVEC*RAD2)/(1.0_dp+WAVEC*RAD1))**2 * DEXP(-2.0_dp*ALPHA*(RAD2-RAD1))
      
      ! Test on acceptance
      if (QUOT < 1) THEN
        call GENRAN(rannumb)
        MCSCHRITT = dble(rannumb) < QUOT
      else
        MCSCHRITT = .TRUE.
      end if
      
      if (MCSCHRITT) THEN
        RE = RNEU
        MCOUNT = MCOUNT + 1
      else
        RNEU = RE
      end if

    end do lrunpre
    
    write(35,*)'STEPMAX = ',STEPMAX
    write(35,*)'prerun: MCPRE= ',MCPRE,' acc. ratio = ', 100.*DBLE(MCOUNT)/DBLE(MCPRE),' %  '
    MCOUNT = 0
    COUNTU = 0
    COUNTD = 0
    
    !
    !  MC loop: main run after thermalizing
    lrun:do IMC=1,MCMAX
        
      do ix=1,3
        ! Shift position at random within +-STEPMAX/2
        call GENRAN(rannumb)
        RDIF = (rannumb-0.5_dp)*STEPMAX
        RNEU(ix) = RE(ix)+RDIF
      end do
        
      ! Calculate wave function ratio psi=(1+c*r)exp(-alpha*r)
      RAD1 = DSQRT(RE(1)**2+RE(2)**2+RE(3)**2)
      RAD2 = DSQRT(RNEU(1)**2+RNEU(2)**2+RNEU(3)**2)
      QUOT=((1.0_dp+WAVEC*RAD2)/(1.0_dp+WAVEC*RAD1))**2* DEXP(-2.0_dp*ALPHA*(RAD2-RAD1))
        
      ! Test on acceptance
      if (QUOT < 1) THEN
        call GENRAN(rannumb)
        MCSCHRITT = dble(rannumb) < QUOT
        if (MCSCHRITT) COUNTU = COUNTU +1
      else
        MCSCHRITT = .TRUE.
        COUNTD = COUNTD + 1
      end if
        
      if (MCSCHRITT) THEN
        RE = RNEU
        MCOUNT = MCOUNT + 1
      else
        RNEU = RE
      end if
        
      RADNEU = DSQRT(RE(1)**2 + RE(2)**2 + RE(3)**2)
        
      ! write (*,*)'RADNEU = ',RADNEU
      ! Update of observables
      if (RADNEU .LT. EMACH) THEN
        LOCEN = -0.5_dp*ALPHA**2 + WAVEC**2 + ALPHA*WAVEC + 3.0_dp*(ALPHA-WAVEC)/2.0_dp/EMACH
        EK = 0.0_dp
        EP = 0.0_dp
      else if (DABS(RADNEU*WAVEC+1) .LT. EMACH) THEN
        EK = -0.5_dp*ALPHA**2
        EK = EK + ALPHA - WAVEC*(1.0_dp+2.0_dp*(ALPHA+WAVEC**2))
        EP = -1.0_dp/RADNEU
        LOCEN = EK + EP
      else
        EK = -0.5_dp*ALPHA**2
        EK = EK + ALPHA/RADNEU-WAVEC*(1.0_dp-ALPHA*RADNEU)/(1.0_dp+WAVEC*RADNEU)/RADNEU
        EP = -1.0_dp/RADNEU
        LOCEN = EK + EP
      end if
        
      !  ERWKIN and ERWPOT miss the correction close to the nucleus
      ERWKIN = dble(IMC-1)/dble(IMC)*ERWKIN +EK/dble(IMC)
      ERWPOT = dble(IMC-1)/dble(IMC)*ERWPOT +EP/dble(IMC)
      call DENSITY1D
      AVRHORAD(1:NRHO)=AVRHORAD(1:NRHO)*dble(IMC-1)/dble(IMC)+RHORAD(1:NRHO)/dble(IMC)
      ERWEN = dble(IMC-1)/dble(IMC)*ERWEN+LOCEN/dble(IMC)
      if (IMC == 1) THEN
        VAREN = 0.0_dp
      else
        VAREN = dble(IMC-1)/dble(IMC)*VAREN + 1/dble(IMC-1)*(ERWEN-LOCEN)**2
      end if
    end do lrun
    
    work1(n2) = VAREN
    work2(n2) = ERWEN
    write(35,35)'main run: MCMAX= ',MCMAX,' acc. ratio = ', 100.*dble(MCOUNT)/dble(MCMAX),' %  '
    write(35,*) 'downhill steps, towards higher probability, COUNTS in %= ', &
               100.D0*COUNTD/dble(MCMAX),' % '
    write(35,*) 'uphill steps, towards lower probability, COUNTS in %= ',100.D0*COUNTU/dble(MCMAX),' % '
    write(35,*) 'ALPHA = ',ALPHA,' WAVEC = ',WAVEC
    write(35,*) 'energy+0.5*ALPHA**2 = ',ERWEN+0.5*ALPHA**2
    write(35,*) 'ERWKIN = ',ERWKIN,'   ERWPOT = ',ERWPOT
    write(35,*) 'ERWEN = ',ERWEN,'  VAREN = ',VAREN
    write(35,*)
    write(35,*)

    write(*,*) 'Done lnxfa2: ', n2, ' from ', NXFA2
  end do lnxfa2

  do n3=1,NXFA2+1
    ! Cut off variance above 0.01 because data plot
    ! gets too large spread in z-values
    if (work1(n3) .GT. 0.01) work1(n3) = 0.05
      write(38,25) ALPHA,WAVEC0 + (n3-1)*DWAVEC,work1(n3)
      write(39,25) ALPHA,WAVEC0 + (n3-1)*DWAVEC,work2(n3)
    end do

    write(*,*) '*** Done lnxfa1: ', n1, ' from ', NXFA1
    write(*,*)
  end do lnxfa1

  do n3=1,NRHO
   write(47,*) n3,AVRHORAD(n3)
  end do

  close(35)
  close(38)
  close(39)

25    format(1x,2f7.3,e12.4)
35    format(1x,a,i11,a,f7.3,a)

end