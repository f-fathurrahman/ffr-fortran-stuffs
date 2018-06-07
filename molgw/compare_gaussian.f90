!=========================================================================
function compare_gaussian(g1,g2) result(same_gaussian)
  use m_gaussian
 implicit none
 logical                   :: same_gaussian
 type(gaussian),intent(in) :: g1,g2
!===== 
!===== 

 same_gaussian = .TRUE.

 if( g1%am /= g2%am                          ) same_gaussian = .FALSE.
 if( g1%nx /= g2%nx                          ) same_gaussian = .FALSE.
 if( g1%ny /= g2%ny                          ) same_gaussian = .FALSE.
 if( g1%nz /= g2%nz                          ) same_gaussian = .FALSE.
 if( ABS(g1%alpha - g2%alpha) > 1.0e-5       ) same_gaussian = .FALSE.
 if( ANY(ABS(g1%x0(:) - g2%x0(:)) > 1.0e-5 ) ) same_gaussian = .FALSE.

end function compare_gaussian
