!=========================================================================
function ank(n,k)
  use m_definitions
  implicit none

  integer,intent(in) :: n,k
  real(dp)           :: ank
  !=====
  integer  :: i
  !=====

  ank   = 1.0_dp
  do i=n,k+1,-1
    ank   = ank   * REAL(i,dp)
  enddo

end function ank
