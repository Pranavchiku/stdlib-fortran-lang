program example_log_factorial
  use stdlib_kinds, only: int64
  use stdlib_specialfunctions_gamma, only: lf => log_factorial
  implicit none
  integer :: n

  n = 10
  print *, lf(n)
  if (abs(lf(n) - 15.1044130) > 1e-8) error stop

! 15.1044130

  print *, lf(35_int64)
  if (abs(lf(35_int64) - 92.1361771) > 1e-8) error stop

! 92.1361771
end program example_log_factorial
