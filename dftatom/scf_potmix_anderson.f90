SUBROUTINE scf_potmix_anderson( Nr, n_orb, x0, max_iter, energy_crit, d, alpha, eps, xout )
  USE m_dft_data, ONLY: dft_data_t
  IMPLICIT NONE 
  ! Finds "x" so that R(x) = 0, uses x0 as the initial estimate
  INTEGER :: Nr
  INTEGER :: n_orb
  REAL(8), INTENT(in) :: x0(Nr)
  INTEGER, INTENT(in) :: max_iter
  LOGICAL, INTENT(in) :: energy_crit
  TYPE(dft_data_t), INTENT(inout) :: d ! Data passed to "F"
  REAL(8), INTENT(in) :: alpha
  REAL(8), INTENT(in) :: eps
  REAL(8) :: xout(Nr)
  REAL(8), PARAMETER :: pi=4.d0*atan(1.d0)
  REAL(8) :: x_i(Nr), x1_i(Nr), R_i(Nr), R1_i(Nr), delta_R(Nr), delta_x(Nr)
  REAL(8) :: beta
  REAL(8) :: sn, sd
  REAL(8) :: ks_energies(n_orb)
  REAL(8) :: x_i_norm, R_i_norm
  REAL(8) :: err_old, err, L2_err
  INTEGER :: i
  REAL(8) :: ss

  x_i = x0
  IF( energy_crit ) THEN 
    ks_energies = d%ks_energies
    err_old = 1d12
  ENDIF 
  
  DO i = 1, max_iter
    !R_i <- R(x_i, i, d)
    CALL ks_potential_map( Nr, x_i, i, d, R_i )
    !
    IF( energy_crit ) THEN 
      ! L2 norm of the "input" potential:
      CALL integrate_trapz_7( Nr, d%Rp, x_i**2 * d%R**2, ss )
      x_i_norm = sqrt( 4*pi*ss )
      ! L2 norm of the "output-input" potential:
      CALL integrate_trapz_7( Nr, d%Rp, R_i**2 * d%R**2, ss )
      R_i_norm = sqrt( 4*pi*ss )
      !
      IF( x_i_norm < 1d-12 ) x_i_norm = 1d-12
      !
      L2_err = R_i_norm / x_i_norm
      err = maxval( abs(ks_energies - d%ks_energies) )
      !
      print *, i, L2_err, err
      ! Do at least 3 iterations
      IF( i >= 3 .and. L2_err < 5d-5) THEN 
        IF( err < eps .and. err_old < eps ) THEN 
          xout = x_i
          RETURN 
        ENDIF 
      ENDIF 
      ks_energies = d%ks_energies
      err_old = err
      !
    ENDIF 

    IF( i > 1 ) THEN 
      delta_x = x_i - x1_i
      delta_R = R_i - R1_i
    ENDIF 
    x1_i = x_i
    R1_i = R_i
    x_i = x_i + alpha * R_i
    IF( i > 1 ) THEN 
      sn = sum(R_i * delta_R * d%R**2)
      sd = sum(delta_R**2 * d%R**2)
      beta = sn / sd
      x_i = x_i - beta * (delta_x + alpha * delta_R)
    ENDIF 
  ENDDO 
  xout = x_i
  IF( energy_crit ) CALL stop_error("SCF didn't converge")
END SUBROUTINE 

