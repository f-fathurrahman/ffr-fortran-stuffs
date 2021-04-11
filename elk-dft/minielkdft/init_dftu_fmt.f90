SUBROUTINE init_dftu_fmt()
  USE modmain
  USE moddftu, ONLY: &
               dmatmt, itmfix, idftu, dftu, tvmatmt, ntmfix, ndftu, lmaxdm, &
               lmmaxdm, ftmtype, tvmmt, vmftm, alphadu, engyadu, vmatmt
  IMPLICIT NONE 
  INTEGER :: i, is, ia, l, ias

  !-------------------------------------------------!
  !     DFT+U and fixed tensor moment variables     !
  !-------------------------------------------------!
  IF( (dftu /= 0) .or. (ftmtype /= 0)) THEN 
    !
    ! density matrix elements in each muffin-tin
    IF( allocated(dmatmt) ) DEALLOCATE( dmatmt )
    ALLOCATE(dmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
    !
    ! potential matrix elements in each muffin-tin
    IF( allocated(vmatmt) ) DEALLOCATE(vmatmt)
    ALLOCATE( vmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot) )
    !
    ! zero the potential matrix
    vmatmt(:,:,:,:,:) = 0.d0
    ! require the potential matrix elements be calculated
    tvmatmt = .true.
    ! flags for non-zero muffin-tin potential matrices
    IF( allocated(tvmmt) ) DEALLOCATE(tvmmt)
    ALLOCATE( tvmmt(0:lmaxdm,natmtot) )
    tvmmt(:,:) = .false.
    ! require second-variational eigenvectors
    tevecsv=.true.
  ENDIF 
  !
  IF( dftu /= 0 ) THEN 
    ! DFT+U energy for each atom
    IF( allocated(engyadu) ) DEALLOCATE( engyadu )
    ALLOCATE( engyadu(natmmax,ndftu) )
    ! interpolation constants (alpha)
    IF( allocated(alphadu) ) DEALLOCATE( alphadu )
    ALLOCATE( alphadu(natmmax,ndftu))
    ! flag the muffin-tin potential matrices which are non-zero
    DO i = 1,ndftu
      is = idftu(1,i)
      IF( is > nspecies) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(init0): invalid species number : ",I8)') is
        WRITE(*,*)
        STOP 
      ENDIF 
      l = idftu(2,i)
      DO ia = 1,natoms(is)
        ias = idxas(ia,is)
        tvmmt(l,ias) = .true.
      ENDDO 
    ENDDO 
  ENDIF 
  !
  IF( ftmtype /= 0 ) THEN 
    ! allocate and zero the fixed tensor moment potential array
    IF( allocated(vmftm) ) DEALLOCATE(vmftm)
    ALLOCATE( vmftm(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot) )
    vmftm(:,:,:,:,:) = 0.d0
    ! flag the muffin-tin potential matrices which are non-zero
    DO i = 1,ntmfix
      is = itmfix(1,i)
      ia = itmfix(2,i)
      ias = idxas(ia,is)
      l = itmfix(3,i)
      tvmmt(l,ias) = .true.
    ENDDO 
  ENDIF 

END SUBROUTINE 

