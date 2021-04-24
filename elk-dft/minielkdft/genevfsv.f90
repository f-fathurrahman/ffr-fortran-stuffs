SUBROUTINE genevfsv()
  USE m_states, ONLY: nstsv, nstfv, evalsv
  USE m_spin, ONLY: nspnfv
  USE m_hamiltonian, ONLY: nmatmax
  USE m_kpoints, ONLY: nkpt
  use modmpi, ONLY: np_mpi, lp_mpi, mpicom, mpi_double_precision, ierror
  USE m_misc, ONLY: filext
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,lp
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: evalfv(:,:)
  COMPLEX(8), ALLOCATABLE :: evecfv(:,:,:),evecsv(:,:)
  ! begin parallel loop over k-points
  ALLOCATE(evalfv(nstfv,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  DO ik=1,nkpt
    ! distribute among MPI processes
    IF(mod(ik-1,np_mpi).ne.lp_mpi) cycle
    ! solve the first- and second-variational eigenvalue equations
    CALL eveqn(ik,evalfv,evecfv,evecsv)
    ! write the eigenvalues/vectors to file
    CALL putevalfv(filext,ik,evalfv)
    CALL putevalsv(filext,ik,evalsv(:,ik))
    CALL putevecfv(filext,ik,evecfv)
    CALL putevecsv(filext,ik,evecsv)
  ENDDO 
  DEALLOCATE(evalfv,evecfv,evecsv)
  ! broadcast eigenvalue array to every MPI process
  IF(np_mpi > 1) THEN 
    DO ik=1,nkpt
      lp=mod(ik-1,np_mpi)
      CALL mpi_bcast(evalsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
    ENDDO 
  ENDIF 
  RETURN 
END SUBROUTINE 

