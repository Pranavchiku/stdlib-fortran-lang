program example_log_gamma
  use stdlib_kinds, only: dp
  use stdlib_specialfunctions_gamma, only: log_gamma
  implicit none

  integer :: i
  real :: x
  real(dp) :: y
  complex :: z
  complex(dp) :: z1

  i = 10
  x = 8.76
  y = x
  z = (5.345, -3.467)
  z1 = z
  print *, log_gamma(i)     !default single precision output
!12.8018274
  if (abs(log_gamma(i) - 12.8018274) > 1.0e-8) error stop

  print *, log_gamma(x)     !intrinsic function call
  if (abs(log_gamma(x) - 10.0942659) > 1.0e-6) error stop

!10.0942659

  print *, log_gamma(y)     !intrinsic function call
  if (abs(log_gamma(y) - 10.094266012080043_dp) > 1.0e-12) error stop

!10.094265528673880

  print *, log_gamma(z)     !same kind as input
  if (abs(log_gamma(z) - (2.56165648, -5.73382425)) > 1.0e-6) error stop

!(2.56165648, -5.73382425)

  print *, log_gamma(z1)
  if (abs(log_gamma(z1) - (2.5616571312593273_dp, -5.7338246618229434_dp)) > 1.0e-12) error stop

!(2.5616575105114614, -5.7338247782852498)
  

end program example_log_gamma
