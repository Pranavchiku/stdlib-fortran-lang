program example_linspace_complex
  use stdlib_math, only: linspace
  use stdlib_kinds, only: dp
  implicit none

  complex(dp) :: start = cmplx(10.0_dp, 5.0_dp, kind=dp)
  complex(dp) :: end = cmplx(-10.0_dp, 15.0_dp, kind=dp)

  complex(dp) :: z(11)

  z = linspace(start, end, 11)
  print *, "z = ", z

  print *, abs(sum(z))
  if (abs(abs(sum(z)) - 1.10000000000000000e+02) > 1e-8) error stop
end program example_linspace_complex
