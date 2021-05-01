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
  INTEGER :: is,ia,ias
  INTEGER :: l,i
  REAL(8) :: ts0,ts1  
  
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

  CALL symmetry() ! crystal symmetry set up

  CALL init_radial_meshes()

  CALL init_charges_states()

  CALL init_gvector_arrays()

  CALL init_atoms_cores()

  CALL init_chgden_pot_xc()

  CALL init_forces()
  
  CALL init_dftu_fmt()
  
  !-----------------------!
  !     miscellaneous     !
  !-----------------------!
  !
  ! determine nuclear radii and volumes
  CALL nuclei()
  !
  ! determine the nuclear-nuclear energy
  CALL energynn()
  !
  ! get smearing function description
  CALL getsdata(stype,sdescr)
  !
  ! get mixing type description
  CALL getmixdata(mixtype,mixdescr)
  !
  ! generate the spherical harmonic transform (SHT) matrices
  CALL genshtmat
  !
  ! allocate 1D plotting arrays
  IF( allocated(dvp1d)) DEALLOCATE(dvp1d)
  ALLOCATE( dvp1d(nvp1d))
  IF(allocated(vplp1d)) DEALLOCATE(vplp1d)
  ALLOCATE( vplp1d(3,npp1d) )
  IF( allocated(dpp1d) ) DEALLOCATE(dpp1d)
  ALLOCATE( dpp1d(npp1d) )
  !
  ! initial self-consistent loop number
  iscl = 1
  tlast = .false.
  !
  ! set the Fermi energy to zero
  efermi = 0.d0
  !
  ! set the temperature from the smearing width
  tempk = swidth/kboltz
  
  CALL timesec(ts1)
  timeinit = timeinit + ts1 - ts0

  RETURN 
END SUBROUTINE 


