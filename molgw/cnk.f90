!=========================================================================
function cnk(n,k)
  use m_definitions
 implicit none

 integer,intent(in) :: n,k
 real(dp)           :: cnk
!=====
 integer  :: i
 real(dp) :: num,denom
!=====

 num   = 1.0_dp
 denom = 1.0_dp
 do i=0,k-1
   num   = num   * REAL(n-i,dp)
   denom = denom * ( i + 1.0_dp)
 enddo
 cnk = num / denom

end function cnk
