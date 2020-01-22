module constant

  implicit none

  integer, parameter :: dp=8
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)
  real(dp), parameter :: gamma=1.4_dp
  real(dp), parameter :: c=2.99792458e8_dp
  real(dp), dimension(1), parameter :: gauss_point1=(/0.0_dp/)
  real(dp), dimension(1), parameter :: gauss_weight1=(/2.0_dp/)
  real(dp), dimension(2), parameter :: gauss_point2=(/-0.5773502691896257_dp,0.5773502691896257_dp/)
  real(dp), dimension(2), parameter :: gauss_weight2=(/1.0_dp,1.0_dp/)
  real(dp), dimension(3), parameter :: gauss_point3=(/0.0_dp,-0.7745966692414834_dp,0.7745966692414834_dp/)
  real(dp), dimension(3), parameter :: gauss_weight3=(/0.8888888888888888_dp,0.5555555555555556_dp,0.5555555555555556_dp/)
  real(dp), dimension(4), parameter :: gauss_point4=(/-0.3399810435848563_dp,0.3399810435848563_dp,-0.8611363115940526_dp, &
       0.8611363115940526_dp/)
  real(dp), dimension(4), parameter :: gauss_weight4=(/0.6521451548625461_dp,0.6521451548625461_dp,0.3478548451374538_dp, &
       0.3478548451374538_dp/)
  real(dp), dimension(5), parameter :: gauss_point5=(/0.0_dp,-0.5384693101056831_dp,0.5384693101056831_dp,-0.9061798459386640_dp, &
       0.9061798459386640_dp/)
  real(dp), dimension(5), parameter :: gauss_weight5=(/0.5688888888888889_dp,0.4786286704993665_dp,0.4786286704993665_dp, &
       0.2369268850561891_dp,0.2369268850561891_dp/)
  real(dp), dimension(6), parameter :: gauss_point6=(/0.6612093864662645_dp,-0.6612093864662645_dp,-0.2386191860831969_dp, &
       0.2386191860831969_dp,-0.9324695142031521_dp,0.9324695142031521_dp/)
  real(dp), dimension(6), parameter :: gauss_weight6=(/0.3607615730481386_dp,0.3607615730481386_dp,0.4679139345726910_dp, &
       0.4679139345726910_dp,0.1713244923791704_dp,0.1713244923791704_dp/)
  real(dp), dimension(6), parameter :: size_stencil=(/1.5_dp,1.5_dp,1.5_dp,1.5_dp,1.5_dp,1.5_dp/)

end module constant
