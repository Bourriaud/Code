module constant

  implicit none

  integer, parameter :: dp=8
  real(dp), parameter :: eps=10.0_dp**(-2.0_dp*dp+1)
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)

end module constant
