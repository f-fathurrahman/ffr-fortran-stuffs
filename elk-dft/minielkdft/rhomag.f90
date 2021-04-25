SUBROUTINE rhomag()
  USE m_atomic, ONLY: natmtot, idxis
  USE m_spin, ONLY: spinpol, nspnfv, ndmag, ncmag
  USE m_mt_rad_am, ONLY: npcmt, nrcmti, nrcmt, npmtmax, lmmaxapw
  USE m_density_pot_xc, ONLY: rhomt, rhoir, magir, magmt
  USE m_apwlo, ONLY: apwordmax
  USE m_states, ONLY: nstfv, nstsv, occsv
  USE m_kpoints, ONLY: nkpt, vkl, wkpt
  USE m_gvectors, ONLY: ngtot
  USE m_gkvectors, ONLY: ngkmax, ngk, igkig, vgkc, vgkl, gkc, sfacgk
  USE m_hamiltonian, ONLY: nmatmax
  USE modmpi, ONLY: np_mpi, mpicom, mpi_double_precision, mpi_in_place, lp_mpi, mpi_sum, ierror
  USE m_misc, ONLY: filext
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,ispn,idm
  INTEGER :: is,ias,n
  ! automatic arrays
  INTEGER(8) lock(natmtot)
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
  
  ! initialise the OpenMP locks
  DO ias=1,natmtot
    CALL omp_init_lock(lock(ias))
  ENDDO 

  ! set the charge density and magnetisation to zero
  DO ias=1,natmtot
    is=idxis(ias)
    rhomt(1:npcmt(is),ias)=0.d0
  ENDDO 
  rhoir(:)=0.d0
  DO idm=1,ndmag
    DO ias=1,natmtot
      is=idxis(ias)
      magmt(1:npcmt(is),ias,idm)=0.d0
    ENDDO 
    magir(:,idm)=0.d0
  ENDDO 

  ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))

  DO ik=1,nkpt
    ! distribute among MPI processes
    IF(mod(ik-1,np_mpi).ne.lp_mpi) cycle
    ! get the eigenvectors from file
    CALL getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    CALL getevecsv(filext,ik,vkl(:,ik),evecsv)
    ! find the matching coefficients
    DO ispn=1,nspnfv
      CALL match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
       sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
    ENDDO 
    ! add to the density and magnetisation
    CALL rhomagk(ngk(:,ik),igkig(:,:,ik),lock,wkpt(ik),occsv(:,ik),apwalm, &
     evecfv,evecsv)
  ENDDO 

  DEALLOCATE(apwalm,evecfv,evecsv)

  ! destroy the OpenMP locks
  DO ias=1,natmtot
    CALL omp_destroy_lock(lock(ias))
  ENDDO 

  ! convert muffin-tin density/magnetisation to spherical harmonics
  CALL rhomagsh()

  ! symmetrise the density
  CALL symrf(nrcmt,nrcmti,npcmt,npmtmax,rhomt,rhoir)

  ! convert the density from a coarse to a fine radial mesh
  CALL rfmtctof(rhomt)

  ! symmetrise the magnetisation
  IF(spinpol) CALL symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,npmtmax,magmt,magir)

  ! convert the magnetisation from a coarse to a fine radial mesh
  DO idm=1,ndmag
    CALL rfmtctof(magmt(:,:,idm))
  ENDDO 

  ! add densities from each process and redistribute
  IF(np_mpi > 1) THEN 
    n=npmtmax*natmtot
    CALL mpi_allreduce(mpi_in_place,rhomt,n,mpi_double_precision,mpi_sum,mpicom, &
     ierror)
    CALL mpi_allreduce(mpi_in_place,rhoir,ngtot,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    IF(spinpol) THEN 
      n=n*ndmag
      CALL mpi_allreduce(mpi_in_place,magmt,n,mpi_double_precision,mpi_sum, &
       mpicom,ierror)
      n=ngtot*ndmag
      CALL mpi_allreduce(mpi_in_place,magir,n,mpi_double_precision,mpi_sum, &
       mpicom,ierror)
    ENDIF 
  ENDIF 

  ! synchronise MPI processes
  CALL mpi_barrier(mpicom,ierror)

  ! add the core density and magnetisation to the total
  CALL rhocore()

  ! calculate the charges
  CALL charge()

  ! normalise the density
  CALL rhonorm()

  ! calculate the moments
  IF(spinpol) CALL moment()

  RETURN 
END SUBROUTINE 

