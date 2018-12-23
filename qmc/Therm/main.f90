PROGRAM THERM
  USE m_random
  USE m_position
  USE m_output, ONLY: FE1MAX,FE2MAX,FE3MAX,PRONAME,OUT3D

  ! NE = number of particles
  IMPLICIT NONE
  INTEGER, PARAMETER :: NXFA1=0,NXFA2=9
  REAL(dp) :: X,Y,Z,STEPMAX,ALPHA,DALPHA,ALPHA0,ALPHA1,BETA
  REAL(dp) :: DLENGTH,LENGTH0,LENGTH1,DUMMY
  REAL(dp) :: RDIF,QUOT
  REAL(dp) :: LOCEN,ERWEN,VAREN,ERWKIN,VARKIN,ERWPOT,VARPOT
  REAL(dp) :: AVRHO(FE3MAX,FE3MAX,FE3MAX)
  LOGICAL :: MCSCHRITT,SWIRHO
  CHARACTER(7) :: FELDNAM
  ! Local variables
  INTEGER :: n1,n2,i,k,ny,nz,M
  REAL(dp) :: RAD1,RAD2,RADNEU,EK,ep,rannumb
  REAL(dp) :: WORK1(NXFA2+1),WORK2(NXFA2+1),WORK
  LOGICAL :: l

  PRONAME = 'thermod'

  ! Controll output
  OPEN(unit=35,file=PRONAME//'_MC10.OUT',position='append', status='unknown')
  
  WRITE(35,*) 'MC7 and MC6 have totally unbiased boundary condition'
  WRITE(35,*) 'totally unbiased is still biased;instead MC8 rejects'
  WRITE(35,*) 'any move outside the volume which is really unbiased'

  flush(35)

  ! Contour plot output  ALPHA and WAVEC for "xfarbe" tool on a NXFA*NXFA grid
  open(unit=38,file=PRONAME//'erw_yukawa10.dat',position='append',status='unknown')
  open(unit=39,file=PRONAME//'var_yukawa10.dat',position='append',status='unknown')

  ! Number of electrons and MC steps
      NE = 100
      IF (NE .GT. NEMAX) THEN
        write (*,*)'NE= ',NE,' larger than NEMAX= ',NEMAX
        STOP
      ENDIF
      MCPRE = 100000
      MCMAX = 1000000
      NDIV = 21    ! NDIV must be odd
  
  call INITRAN()

!  Start data:
      lparam1:do m=1,3
      lparam2:do n1 = 1,NXFA1+1
      lparam3:do n2 = 1,NXFA2+1
!  Parameter ALPHA of inverse temperature and edge length LENGTH
      SWIRHO = .true.
      ALPHA0=0.50D0
      ALPHA1=10.0D0
      LENGTH0=+1.00D0
      LENGTH1=+50.0D0
      DALPHA=(ALPHA1-ALPHA0)/DBLE(NXFA1+1)
      DLENGTH=(LENGTH1-LENGTH0)/DBLE(NXFA2+1)
      ALPHA=ALPHA0+(n1-1)*DALPHA
!  adjust LENGTH to logarithmic plot
      LENGTH=(0.01+(n2-1)*0.01)*10**m
!      LENGTH=LENGTH0+(n2-1)*DLENGTH
!  Maximum step width, KONTUZ: Check with acceptance ratio!
      STEPMAX = 9.0D-2*LENGTH
!  Random initial electron positions
      do k=1,NE
        do i=1,3
          call GENRAN(rannumb)
          RE(i,k) =  LENGTH*(rannumb-0.5)
          RNEU(i,k) = RE(i,k)
        end do
!      write(*,*)'k= ',k,(RE(i,k),i=1,3)
      end do
      do i=1,NE
       VNEW(i) = 0.D0
       VDIF(i) = 0.D0
       do k=1,NE
        IF (k .eq. i) cycle
        DISTNEU(1:3,i,k) = RNEU(1:3,i)-RNEU(1:3,k)
        DIST(1:3,i,k) = DISTNEU(1:3,i,k)
        DISTNEU(4,i,k) = dsqrt(sum (DISTNEU(1:3,i,k)**2))
        DIST(4,i,k) = DISTNEU(4,i,k)
        VNEW(i) = VNEW(i) + 1.D0/DISTNEU(4,i,k)*dexp(-DISTNEU(4,i,k))
       end do
      end do
!  Counts the acceptance number
      MCOUNT = 0
!  Observables
      RHO(1:NDIV,1:NDIV,1:NDIV) = 0.D0
      AVRHO(1:NDIV,1:NDIV,1:NDIV) = 0.D0
      LOCEN = 0.D0
      ERWEN = 0.D0
      VAREN = 0.D0
      ERWKIN = 0.D0
      VARKIN = 0.D0
      ERWPOT = 0.D0
      VARPOT = 0.D0
      l = .true.
!
!
!  MC loop: prerun for thermalizing,KONTUZ: does not change the bad
!           sampling of energy!!!
      lrunpre:do IMC=2,MCPRE
       lelrunpre:do IE=1,NE
        do i=1,3
!  Shift position at random within +-STEPMAX/2
          call GENRAN(rannumb)
          RDIF = (rannumb-0.5)*STEPMAX
          RNEU(i,IE) = RE(i,IE)+RDIF
!CC  reflect particle at wall if it crosses
!C          if (RNEU(i,IE) > 0.5D0*LENGTH)
!C     )              RNEU(i,IE)=LENGTH-RNEU(i,IE)
!C          IF (RNEU(i,IE) < -0.5D0*LENGTH)
!C     )              RNEU(i,IE)=-LENGTH-RNEU(i,IE)
!CC  scatter particle at a rough wall
!C          if (DABS(RNEU(i,IE)) > 0.5D0*LENGTH) THEN
!C           call GENRAN(rannumb)
!C            if (RNEU(i,IE) > 0.5D0*LENGTH)
!C     &        RNEU(i,IE)=LENGTH/2.D0-ABS(rannumb-0.5)*STEPMAX
!C     &        RNEU(i,IE)=LENGTH*(1/2.D0-rannumb) ! totally unbiased
!C            if (RNEU(i,IE) < -0.5D0*LENGTH)
!C     &        RNEU(i,IE)=-LENGTH/2.D0+ABS(rannumb-0.5)*STEPMAX
!C     &        RNEU(i,IE)=-LENGTH*(1/2.D0-rannumb) ! totally unbiased
!   Apply strict cut-off boundary conditions
          if (DABS(RNEU(i,IE)) > 0.5D0*LENGTH) THEN
           MCSCHRITT = .false.
           RNEU(1:3,IE) = RE(1:3,IE)
           MCOUNT = MCOUNT + 1
           cycle lelrunpre
          end if
        end do
!  Calculate Boltzmann ratio with energy 0.5*sum_k 1/r_ik exp(-r_ik)
!    without term k=i
         call JASEXP(VNEW(IE),VDIF(IE))
         QUOT = DEXP(-ALPHA*VDIF(IE))
!  Test on acceptance
         if (QUOT < 1) THEN
!           MCSCHRITT = DBLE(ZBQLU01(DUMMY)) .LT. QUOT
           call GENRAN(rannumb)
           MCSCHRITT = DBLE(rannumb) < QUOT
         else
           MCSCHRITT = .true.
         end if
         if (MCSCHRITT) then
             RE(1:3,IE) = RNEU(1:3,IE)
             MCOUNT = MCOUNT + 1
         else
             RNEU(1:3,IE) = RE(1:3,IE)
         end if
       end do lelrunpre
      end do lrunpre
      WRITE(35,*) 'STEPMAX = ',STEPMAX
      WRITE(35,*) 'prerun: MCPRE= ',MCPRE,' acc. ratio = ', 100.*DBLE(MCOUNT)/DBLE(NE*MCPRE),' %  '
  flush(35)
MCOUNT = 0
!
!
!
!  MC loop: main run after thermalizing
      lrun:do IMC=2,MCMAX
       lelrun:do IE=1,NE
        do i=1,3
!  Shift position at random within +-STEPMAX/2
          call GENRAN(rannumb)
          RDIF = (rannumb-0.5)*STEPMAX
          RNEU(i,IE) = RE(i,IE)+RDIF
!CC  Reflect particle at wall if it crosses
!C          if (RNEU(i,IE) > 0.5D0*LENGTH)
!C     )              RNEU(i,IE)=LENGTH-RNEU(i,IE)
!C          if (RNEU(i,IE) < -0.5D0*LENGTH)
!C     )              RNEU(i,IE)=-LENGTH-RNEU(i,IE)
!CC  Scatter particle at a rough wall
!C          if (DABS(RNEU(i,IE)) .gt. 0.5D0*LENGTH) then
!C            call GENRAN(rannumb)
!C            if (RNEU(i,IE) > 0.5D0*LENGTH)
!C     &        RNEU(i,IE)=LENGTH/2.D0-ABS(ZBQLU01(DUMMY)-0.5)*STEPMAX
!C            if (RNEU(i,IE) < -0.5D0*LENGTH)
!C     &        RNEU(i,IE)=-LENGTH/2.D0+ABS(ZBQLU01(DUMMY)-0.5)*STEPMAX
!C     &        RNEU(i,IE)=LENGTH/2.D0-ABS(rannumb-0.5)*STEPMAX
!C     &        RNEU(i,IE)=LENGTH*(1/2.D0-rannumb) ! totally unbiased
!C            if (RNEU(i,IE) .lt. -0.5D0*LENGTH)
!C     &        RNEU(i,IE)=-LENGTH/2.D0+ABS(rannumb-0.5)*STEPMAX
!C     &        RNEU(i,IE)=-LENGTH*(1/2.D0-rannumb) ! totally unbiased
!   Apply strict cut-off boundary conditions
          if (DABS(RNEU(i,IE)) > 0.5D0*LENGTH) THEN
            MCSCHRITT = .false.
            l = .false.
            exit
          end if
        end do
        if (l) then
! Calculate Boltzmann ratio with energy 0.5*sum_k 1/r_ik exp(-r_ik) without term k=i
      call JASEXP(VNEW(IE),VDIF(IE))
         QUOT = DEXP(-ALPHA*VDIF(IE))
!  Test on acceptance
         if (QUOT < 1) then
           call GENRAN(rannumb)
           MCSCHRITT = DBLE(rannumb) < QUOT
         else
           MCSCHRITT = .true.
         end if
        end if
        l = .true.
!  Update of observables.
        if (MCSCHRITT) then
          RE(1:3,IE) = RNEU(1:3,IE)
          ep = VNEW(IE)
          MCOUNT = MCOUNT + 1
        else
          RNEU(1:3,IE) =  RE(1:3,IE)
          ep = VNEW(IE) - VDIF(IE)
        end if
!  E=0.5 sum_ik v_ik, sum i omitted because
!  run averages over particles, so energy per particle is calculated
        LOCEN = LOCEN + 0.5D0*ep
       end do lelrun
!  energy per particle
       LOCEN = LOCEN/DBLE(NE)
       if (IMC == 2) then  !set start values
        ERWEN = LOCEN
        VAREN = 0.D0
       end if
       ERWEN = (IMC-1)/DBLE(IMC)*ERWEN+LOCEN/DBLE(IMC)
       VAREN = (IMC-1)/DBLE(IMC)*VAREN + 1/DBLE(IMC-1)*(ERWEN-LOCEN)**2
       LOCEN = 0.D0
       if (SWIRHO) then
!  density
        call DENSITY(RHO)
        do nz=1,NDIV
         do ny=1,NDIV
          AVRHO(1:NDIV,ny,nz)=(IMC-1)/DBLE(IMC)*AVRHO(1:NDIV,ny,nz) &
                        +  RHO(1:NDIV,ny,nz)/DBLE(IMC)
         end do
        end do
       end if
      end do lrun
!  end MC loop
!
      WORK1(n2) = ERWEN
      WORK2(n2) = VAREN
  write(35,35) 'main run: MCMAX= ', MCMAX, &
               ' acc. ratio = ', 100.*DBLE(MCOUNT)/DBLE(NE*MCMAX),' %  '
      write(35,*)'ALPHA = ',ALPHA,' LENGTH = ',LENGTH
      write(35,*)'energy+0.5*ALPHA**2 = ',ERWEN+0.5*ALPHA**2
      write(35,*)'ERWEN = ',ERWEN,'  VAREN = ',VAREN
      write(35,*)
      write(35,*)
      if (SWIRHO) then
        write (36,*) 'ALPHA = ',ALPHA,' LENGTH = ',LENGTH
        call OUT3D(NDIV,NDIV,NDIV,AVRHO)
      end if
      write(38,45) LENGTH,ERWEN
      write(39,45) LENGTH,VAREN
      end do lparam3
      write(38,fmt="(t3)")
      write(39,fmt="(t3)")
      end do lparam2
      end do lparam1
      CLOSE(35)
      CLOSE(38)
      CLOSE(39)
35    FORMAT(1X,A,I8,A,F7.3,A)
45    FORMAT(t3,F7.3,E12.3)
end program

