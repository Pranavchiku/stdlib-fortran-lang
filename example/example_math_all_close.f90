program example_math_all_close

  use stdlib_math, only: all_close
  real    :: y, NAN
  complex :: z(4, 4), z_(4,4)

  y = -3
  NAN = sqrt(y)
  z = cmplx(1.0, 1.0)
  z_ = cmplx(1.0, 1.0) + cmplx(1.0e-11, 1.0e-11)

  print *, all_close(z_, z)     ! T
  print *, all_close([NAN], [NAN]), all_close([NAN], [NAN], equal_nan=.true.)
! NAN, F, T

end program example_math_all_close
