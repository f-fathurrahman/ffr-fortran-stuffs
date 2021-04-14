
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

SUBROUTINE genapwlofr()
  IMPLICIT NONE 
  
  ! generate the APW radial functions
  CALL genapwfr()
  
  ! generate the local-orbital radial functions
  CALL genlofr()
  
  ! compute the overlap radial integrals
  CALL olprad()
  
  ! compute the Hamiltonian radial integrals
  CALL hmlrad()
  
  RETURN 
END SUBROUTINE 

