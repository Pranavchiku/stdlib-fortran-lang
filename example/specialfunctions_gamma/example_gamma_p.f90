program example_gamma_p
  use stdlib_specialfunctions_gamma, only: rgp => regularized_gamma_p
  implicit none

  print *, rgp(3.0, 5.0)
  if (abs(rgp(3.0, 5.0) - 0.875347972) > 1e-9) error stop

! 0.875347972
end program example_gamma_p
