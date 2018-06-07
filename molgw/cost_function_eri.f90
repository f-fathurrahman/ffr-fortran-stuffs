!=========================================================================
! Rough evaluation of the CPU time to get an ERI as a function of the 
! angular momentum
! 
!=========================================================================
function cost_function_eri(am)
  USE m_definitions, ONLY : dp
 implicit none
 integer,intent(in)  :: am
 real(dp)            :: cost_function_eri
!=====

 cost_function_eri = am**2 + 4.6_dp 

end function cost_function_eri
