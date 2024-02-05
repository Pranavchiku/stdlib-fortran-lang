program example_logspace_rstart_cbase
  use stdlib_math, only: logspace
  use stdlib_kinds, only: dp
  implicit none

  real(dp) :: start = 0.0_dp
  real(dp) :: end = 3.0_dp
  integer, parameter :: n = 4
  complex(dp) :: base = (0.0_dp, 1.0_dp)

  complex(dp) :: z(n) ! complex values raised to real powers result in complex values

  z = logspace(start, end, n, base)
  print *, "z = ", z
  print *, abs(sum(z))
  if (abs(abs(sum(z)) - 1.83697019872102969e-16_dp) > 1e-12) error stop
end program example_logspace_rstart_cbase
