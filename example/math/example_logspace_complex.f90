program example_logspace_complex
  use stdlib_math, only: logspace
  use stdlib_kinds, only: dp
  implicit none

  complex(dp) :: start = (10.0_dp, 5.0_dp)
  complex(dp) :: end = (-10.0_dp, 15.0_dp)

  complex(dp) :: z(11) ! Complex values raised to complex powers results in complex values

  z = logspace(start, end, 11)
  print *, "z = ", z
  print *, abs(sum(z))

  if (abs(abs(sum(z)) - 9.93335211079602432e+09_dp) > 1e-12) error stop
end program example_logspace_complex
