!=========================================================================
pure function smooth_step(mu)
  use m_definitions
  implicit none
  real(dp) :: smooth_step
  real(dp),intent(in) :: mu

  smooth_step = 0.5_dp * mu * ( 3.0_dp - mu**2 )

end function smooth_step
