SUBROUTINE init_dftu_k()

  USE modmain, ONLY: natmtot, nrmtmax
  USE moddftu, ONLY: dftu, fdue, fdufr, lmaxdm
  IMPLICIT NONE 

  !-------------------------!
  !     DFT+U variables     !
  !-------------------------!
  IF( dftu /= 0 ) THEN 
    ! allocate energy arrays to calculate Slater integrals with Yukawa potential
    IF( allocated(fdue) ) DEALLOCATE(fdue)
    ALLOCATE( fdue(0:lmaxdm, natmtot) )
    ! allocate radial functions to calculate Slater integrals with Yukawa potential
    IF( allocated(fdufr) ) DEALLOCATE(fdufr)
    ALLOCATE( fdufr(nrmtmax, 0:lmaxdm, natmtot) )
  ENDIF

END SUBROUTINE 

