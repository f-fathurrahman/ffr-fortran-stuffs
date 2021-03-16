SUBROUTINE calc_E_nl(c, n, l, Z, relat, E_nl)
  IMPLICIT NONE 
  ! Calculates exact energy for the radial Schroedinger/Dirac equations
  REAL(8), INTENT(in) :: c ! speed of light in atomic units
  INTEGER, INTENT(in) :: n, l, Z, relat
  REAL(8) :: E_nl
  ! quantum numbers (n, l), atomic number (z)
  ! relat == 0 ... Schroedinger equation
  ! relat == 2 ... Dirac equation, spin up
  ! relat == 3 ... Dirac equation, spin down
  
  INTEGER :: kappa
  REAL(8) :: beta
  
  IF( .not. (l >= 0) ) CALL stop_error("'l' must be positive or zero")
  IF( .not. (n > l) ) CALL stop_error("'n' must be greater than 'l'")
  IF( l == 0 .and. relat == 3) CALL stop_error("Spin must be up for l==0.")
  IF( relat == 0 ) THEN
    E_nl = - Z**2 / (2.d0 * n**2)
  ELSE 
    IF( relat == 2 ) THEN 
      kappa = -l - 1
    ELSE 
      kappa = l
    ENDIF 
    beta = sqrt(kappa**2 - (Z/c)**2)
    E_nl = c**2/sqrt(1 + (Z/c)**2/(n - abs(kappa) + beta)**2) - c**2
  ENDIF 
  RETURN 
END SUBROUTINE 

