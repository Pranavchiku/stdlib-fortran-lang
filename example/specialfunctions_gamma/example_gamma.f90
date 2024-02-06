program example_gamma
  use stdlib_kinds, only: dp, int64
  use stdlib_specialfunctions_gamma, only: gamma
  implicit none

  integer :: i
  integer(int64) :: n
  real :: x
  real(dp) :: y
  complex :: z
  complex(dp) :: z1

  i = 10
  n = 15_int64
  x = 2.5
  y = 4.3_dp
  z = (2.3, 0.6)
  z1 = (-4.2_dp, 3.1_dp)

  print *, gamma(i)              !integer gives exact result
! 362880
  if (gamma(i) /= 362880) error stop

  print *, gamma(n)
! 87178291200
  if (gamma(n) /= 87178291200_8) error stop

  print *, gamma(x)              ! intrinsic function call
! 1.32934034
  if (abs(gamma(x) - 1.32934034) > 1.0E-6) error stop

  print *, gamma(y)              ! intrinsic function call
! 8.8553433604540341
  if (abs(gamma(y) - 8.8553433604540341_dp) > 1.0E-12) error stop

  print *, gamma(z)
! (0.988054395, 0.383354813)
  if (abs(gamma(z) - (0.988054395, 0.383354813)) > 1.0E-8) error stop

  print *, gamma(z1)
! (-2.78916032990983999E-005, 9.83164600163221218E-006)
  if (abs(gamma(z1) - (-2.78916032990983999E-005_dp, 9.83164600163221218E-006_dp)) > 1.0E-8) error stop
end program example_gamma
