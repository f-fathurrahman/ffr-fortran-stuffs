SUBROUTINE gndstate
  USE modmain
  USE moddftu
  USE modulr
  USE modmpi
  USE modomp

  IMPLICIT NONE 
  ! local variables
  LOGICAL :: exist
  INTEGER :: ik,nwork,n
  REAL(8) :: dv,etp,de,timetot
  ! allocatable arrays
  REAL(8), allocatable :: v(:),work(:)

  ! initialise global variables
  CALL init0()
  CALL init1()

  ! initialise q-vector-dependent variables if required
  IF( xctype(1) < 0 ) THEN 
    WRITE(*,*) '*** Should call init2'
    CALL init2()
  ELSE 
    WRITE(*,*)
    WRITE(*,*) '*** Not calling init2'
  ENDIF

  ! apply strain to the G, k, G+k and q-vectors if required
  ! FIXME: Check this again
  !IF( (istrain < 1) .or. ( istrain > nstrain) ) THEN 
  !  WRITE(*,*)
  !  WRITE(*,*) '*** Not calling straingkq'
  !ENDIF 
  !call straingkq()

  !
  IF( task == 0) trdstate=.false.
  IF( task == 1) trdstate=.true.

  ! only the MPI master process should write files
  ! !!!! REMOVED !!!!
  
  iscl = 0

  IF( trdstate ) THEN 
    STOP 'Reading from file is deactivated'
    ! read the Kohn-Sham potential and fields from file
    !call readstate
    !if (mp_mpi) then
    !  write(6,'("Potential read in from STATE.OUT")')
    !end if
    !if (autolinengy) call readfermi
  ELSE 
    ! initialise the density and magnetisation from atomic data
    CALL rhoinit()
    CALL maginit()
    ! generate the Kohn-Sham potential and magnetic field
    CALL potks(.true.)
    IF( mp_mpi ) THEN 
      WRITE(6,'("Kohn-Sham potential initialised from atomic data")')
    ENDIF 
  ENDIF 

  IF( mp_mpi ) flush(6)

  ! Fourier transform of KS effective potential
  CALL genvsig() 

  ! size of mixing vector
  n = npmtmax*natmtot + ngtot
  IF( spinpol ) n = n + (npcmtmax*natmtot+ngtot)*ndmag
  IF( tvmatmt ) n = n + 2*((lmmaxdm*nspinor)**2)*natmtot

  ! allocate mixing array
  ALLOCATE( v(n) )

  ! determine the size of the mixer work array
  nwork=-1
  CALL mixerifc(mixtype,n,v,dv,nwork,v)
  ALLOCATE( work(nwork) )

  ! initialise the mixer
  iscl = 0
  CALL mixpack(.true.,n,v)
  CALL mixerifc(mixtype,n,v,dv,nwork,work)

  ! set the stop signal to .false.
  tstop = .false.

  ! set last self-consistent loop flag
  tlast = .false.
  etp = 0.d0

  ! begin the self-consistent loop
  IF( mp_mpi ) THEN 
    WRITE(6,*)
    WRITE(6,'("+------------------------------+")')
    WRITE(6,'("| Self-consistent loop started |")')
    WRITE(6,'("+------------------------------+")')
  ENDIF

  DO iscl = 1,maxscl
  
    IF( mp_mpi ) THEN 
      WRITE(6,*)
      WRITE(6,'("+--------------------+")')
      WRITE(6,'("| Loop number : ",I4," |")') iscl
      WRITE(6,'("+--------------------+")')
    ENDIF 
    IF( iscl >= maxscl ) THEN 
      IF( mp_mpi ) THEN 
        WRITE(6,*)
        WRITE(6,'("Reached self-consistent loops maximum")')
      ENDIF 
      WRITE(*,*)
      WRITE(*,'("Warning(gndstate): failed to reach self-consistency after ",I4," loops")') iscl
      tlast = .true.
    ENDIF 
    IF( mp_mpi ) flush(6)
  
    ! reset the OpenMP thread variables
    CALL omp_reset()

    ! generate the core wavefunctions and densities
    CALL gencore()
  
    ! find the new linearisation energies
    CALL linengy()
  
    ! write out the linearisation energies
    IF( mp_mpi ) CALL writelinen()
  
    ! generate the APW and local-orbital radial functions and integrals
    CALL genapwlofr()
  
    ! generate the spin-orbit coupling radial functions
    CALL gensocfr()
  
    ! generate the first- and second-variational eigenvectors and eigenvalues
    CALL genevfsv()
  
    ! find the occupation numbers and Fermi energy
    CALL occupy()
  
    IF( mp_mpi .and. autoswidth ) THEN 
      WRITE(6,*)
      WRITE(6,'("New smearing width : ",G18.10)') swidth
    ENDIF 
  
    IF( mp_mpi ) THEN 
      ! write the occupation numbers to file
      DO ik=1,nkpt
        CALL putoccsv(filext,ik,occsv(:,ik))
      ENDDO 
      ! write eigenvalues to file
      CALL writeeval()
      ! write the Fermi energy to file
      CALL writefermi()
    ENDIF
  
    ! synchronise MPI processes
    CALL mpi_barrier(mpicom,ierror)

    ! generate the density and magnetisation
    CALL rhomag()

    ! DFT+U or fixed tensor moment calculation
    ! .... SKIPPED ....
    ! Write DFT+U
    ! .... SKIPPED ....

    ! compute the Kohn-Sham potentials and magnetic fields
    CALL potks(.true.)

    IF( mp_mpi ) THEN 
      IF( (xcgrad == 3) .and. (c_tb09 /= 0.d0)) THEN 
        WRITE(6,*)
        WRITE(6,'("Tran-Blaha ''09 constant c : ",G18.10)') c_tb09
      ENDIF 
    ENDIF 

    ! pack interstitial and muffin-tin potential and field into one array
    CALL mixpack(.true.,n,v)

    ! mix in the old potential and field with the new
    CALL mixerifc(mixtype, n, v, dv, nwork, work)

    ! make sure every MPI process has a numerically identical potential
    IF(np_mpi > 1) THEN 
      CALL mpi_bcast(v,n,mpi_double_precision,0,mpicom,ierror)
    ENDIF 

    ! unpack potential and field
    CALL mixpack(.false.,n,v)

    ! calculate and add the fixed spin moment effective field (after mixing)
    ! .... SKIPPED ....

    ! Fourier transform Kohn-Sham potential to G-space
    CALL genvsig()

    ! reduce the external magnetic fields if required
    ! .... SKIPPED .... 

    ! compute the energy components
    CALL energy()

    IF( mp_mpi ) THEN 
      ! output energy components
      CALL writeengy(6)
      WRITE(6,*)
      WRITE(6,'("Density of states at Fermi energy : ",G18.10)') fermidos
      WRITE(6,'(" (states/Hartree/unit cell)")')
      WRITE(6,*)
      WRITE(6,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
      WRITE(6,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
      WRITE(6,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
      WRITE(6,'(" at k-point ",I6)') ikgap(3)

      ! write total energy to TOTENERGY.OUT
      WRITE(61,'(G22.12)') engytot
      flush(61)

      ! write DOS at Fermi energy to FERMIDOS.OUT
      WRITE(62,'(G18.10)') fermidos
      flush(62)
    
      ! output charges and moments
      CALL writechg(6)
    
      IF( spinpol ) THEN 
        CALL writemom(6)
        ! write total moment to MOMENT.OUT
        WRITE(63,'(3G18.10)') momtot(1:ndmag)
        flush(63)
        ! write total moment magnitude to MOMENTM.OUT
        WRITE(68,'(G18.10)') momtotm
        flush(68)
      ENDIF 
    
      ! write estimated Kohn-Sham indirect band gap
      WRITE(64,'(G22.12)') bandgap(1)
      flush(64)
    
      ! output effective fields for fixed spin moment calculations
      IF( fsmtype /= 0 ) CALL writefsm(6)
    
      ! check for WRITE file
      inquire(file='WRITE',exist=exist)
      IF( exist ) THEN 
        WRITE(6,*)
        WRITE(6,'("WRITE file exists - writing STATE.OUT")')
        CALL writestate
        OPEN(50,file='WRITE')
        CLOSE(50,status='DELETE')
      ENDIF 
      ! write STATE.OUT file if required
      IF( nwrite >= 1) THEN 
        IF(mod(iscl,nwrite) == 0) THEN 
          call writestate()
          WRITE(6,*)
          WRITE(6,'("Wrote STATE.OUT")')
        ENDIF 
      ENDIF 

      ! write OEP residual
      IF( xctype(1) < 0) THEN 
        WRITE(6,*)
        WRITE(6,'("Magnitude of OEP residual : ",G18.10)') resoep
        WRITE(69,'(G18.10)') resoep
        flush(69)
      ENDIF 
    ENDIF 

    ! exit self-consistent loop if required
    IF( tlast ) GOTO 10

    ! check for convergence
    IF( iscl >= 2) THEN 
      IF( mp_mpi ) THEN 
        WRITE(6,*)
        WRITE(6,'("RMS change in Kohn-Sham potential (target) : ",G18.10," (",G18.10,")")') dv,epspot
        WRITE(65,'(G18.10)') dv
        flush(65)
      ENDIF
      de = abs(engytot - etp)
      IF( mp_mpi ) THEN 
        WRITE(6,'("Absolute change in total energy (target)   : ",G18.10," (",G18.10,")")') de,epsengy
        WRITE(66,'(G18.10)') de
        flush(66)
      ENDIF 
      !
      IF( (dv < epspot) .and. (de < epsengy) ) THEN 
        IF( mp_mpi ) THEN 
          WRITE(6,*)
          WRITE(6,'("Convergence targets achieved")')
        ENDIF 
        tlast=.true.
      ENDIF 
    ENDIF ! check for convergence

    ! average the current and previous total energies and store
    IF (iscl.gt.1) THEN 
      etp = 0.75d0*engytot + 0.25d0*etp
    ELSE 
      etp = engytot
    ENDIF 
  
    ! check for STOP file (only master process)
    IF( mp_mpi ) THEN 
      inquire(file='STOP',exist=exist)
      IF( exist ) THEN 
        WRITE(6,*)
        WRITE(6,'("STOP file exists - stopping self-consistent loop")')
        OPEN(50,file='STOP')
        CLOSE(50,status='DELETE')
        tstop=.true.
        tlast=.true.
      ENDIF 
    ENDIF 

    ! broadcast tlast and tstop from master process to all other processes
    CALL mpi_bcast(tlast,1,mpi_logical,0,mpicom,ierror)
    CALL mpi_bcast(tstop,1,mpi_logical,0,mpicom,ierror)

    ! output the current total CPU time
    timetot = timeinit + timemat + timefv + timesv + timerho + timepot + timefor
    IF( mp_mpi ) THEN 
      WRITE(6,*)
      WRITE(6,'("Time (CPU seconds) : ",F12.2)') timetot
    ENDIF 
  
  ! end the self-consistent loop
  ENDDO 

  10 CONTINUE

  ! synchronise MPI processes
  CALL mpi_barrier(mpicom,ierror)
  IF( mp_mpi ) THEN 
    WRITE(6,*)
    WRITE(6,'("+------------------------------+")')
    WRITE(6,'("| Self-consistent loop stopped |")')
    WRITE(6,'("+------------------------------+")')
    ! write density and potentials to file only if maxscl > 1
    IF( maxscl > 1) THEN 
      CALL writestate()
      WRITE(6,*)
      WRITE(6,'("Wrote STATE.OUT")')
    ENDIF 
  ENDIF 

  ! compute forces if required
  !if (tforce) then
  !  call force()
  !  ! output forces to INFO.OUT
  !  if (mp_mpi) call writeforces(6)
  !end IF

  ! compute the current density and total current if required
  !if (tcden) then
  !  call curden(afieldc)
  !  if (mp_mpi) then
  !    write(6,*)
  !    write(6,'("Total current per unit cell")')
  !    write(6,'(3G18.10)') curtot
  !    write(6,'(" magnitude : ",G18.10)') curtotm
  !  end if
  !end IF

  ! total time used
  timetot = timeinit + timemat + timefv + timesv + timerho + timepot + timefor

  IF( mp_mpi ) THEN 
    ! output timing information
    WRITE(6,*)
    WRITE(6,'("Timings (CPU seconds) :")')
    WRITE(6,'(" initialisation",T40,": ",F12.2)') timeinit
    write(6,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
    WRITE(6,'(" first-variational eigenvalue equation",T40,": ",F12.2)') timefv
    IF( tevecsv ) THEN 
      WRITE(6,'(" second-variational calculation",T40,": ",F12.2)') timesv
    ENDIF 
    WRITE(6,'(" charge density calculation",T40,": ",F12.2)') timerho
    WRITE(6,'(" potential calculation",T40,": ",F12.2)') timepot
    IF( tforce ) THEN 
      WRITE(6,'(" force calculation",T40,": ",F12.2)') timefor
    ENDIF 
    WRITE(6,'(" total",T40,": ",F12.2)') timetot
    WRITE(6,*)
    WRITE(6,'("+----------------------------+")')
    WRITE(6,'("| Elk version ",I1.1,".",I1.1,".",I2.2," stopped |")') version
    WRITE(6,'("+----------------------------+")')

    ! close the INFO.OUT file ! ffr
    !  close(60)

    ! close the TOTENERGY.OUT file
    CLOSE(61)

    ! close the FERMIDOS.OUT file
    CLOSE(62)
    ! close the MOMENT.OUT and MOMENTM.OUT files
    IF( spinpol ) THEN 
      CLOSE(63); close(68)
    ENDIF 
    ! close the GAP.OUT file
    CLOSE(64)
    ! close the RMSDVS.OUT file
    CLOSE(65)
    ! close the DTOTENERGY.OUT file
    CLOSE(66)
    ! close TMDFTU.OUT file
    IF(tmwrite) CLOSE(67)
    ! close the RESIDUAL.OUT file
    IF( xctype(1) < 0 ) CLOSE(69)
  ENDIF 

  deallocate(v,work)

! synchronise MPI processes
call mpi_barrier(mpicom,ierror)

return

end subroutine
!EOC

