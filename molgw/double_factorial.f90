!==========================================
pure function double_factorial(intin)
  use m_definitions
  implicit none
  integer,intent(in) :: intin
  real(dp) :: double_factorial
  !=====
  ! just hard coded for some small integers

 select case (intin)
 case(-1) 
   double_factorial = 1.0_dp
 case( 0) 
   double_factorial = 1.0_dp
 case( 1) 
   double_factorial = 1.0_dp
 case( 2) 
   double_factorial = 2.0_dp
 case( 3) 
   double_factorial = 3.0_dp
 case( 4) 
   double_factorial = 8.0_dp
 case( 5) 
   double_factorial = 15.0_dp
 case( 6) 
   double_factorial = 48.0_dp
 case( 7) 
   double_factorial = 105.0_dp
 case( 8) 
   double_factorial = 384.0_dp
 case( 9) 
   double_factorial = 945.0_dp
 case(10) 
   double_factorial = 3840.0_dp
 case(11) 
   double_factorial = 10395.0_dp
 case(12) 
   double_factorial = 46080.0_dp
 case(13) 
   double_factorial = 135135.0_dp
 case(14) 
   double_factorial = 645120.0_dp
 case(15)
   double_factorial = 2027025.0_dp
 case(16) 
   double_factorial = 10321920.0_dp
 case(17) 
   double_factorial = 34459425.0_dp
 case(18) 
   double_factorial = 185794560.0_dp
 case(19) 
   double_factorial = 654729075.0_dp
 case(20) 
   double_factorial = 3715891200.0_dp
 case(21) 
   double_factorial = 13749310575.0_dp
 case(22)
   double_factorial = 81749606400.0_dp
 case(23) 
   double_factorial = 316234143225.0_dp
 case(25) 
   double_factorial = 7905853580625.0_dp
 case(27) 
   double_factorial = 213458046676875.0_dp
 case(29) 
   double_factorial = 6190283353629375.0_dp
 case(31) 
   double_factorial = 191898783962510625.0_dp
 end select

end function double_factorial