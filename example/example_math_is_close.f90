program example_math_is_close

  use stdlib_math, only: is_close
  real :: x(2) = [1, 2], y, NAN
  logical :: result(2)

  y = -3
  NAN = sqrt(y)

  print *, is_close(x, [real :: 1, 2.1])       ! [T, F]
  result = is_close(x, [real :: 1, 2.1])
  if (.not. result(1)) error stop
  if (result(2)) error stop

  print *, is_close(2.0, 2.1, abs_tol=0.1)    ! T
  if (.not. is_close(2.0, 2.1, abs_tol=0.1)) error stop
  print *, NAN, is_close(2.0, NAN), is_close(2.0, NAN, equal_nan=.true.)   ! NAN, F, F

  if (is_close(NAN, NAN)) error stop
  ! if (is_close(NAN, NAN, equal_nan=.true.)) error stop ! TODO: fix this

  print *, is_close(NAN, NAN), is_close(NAN, NAN, equal_nan=.true.)        ! F, T
  if (is_close(NAN, NAN)) error stop
  if (.not. is_close(NAN, NAN, equal_nan=.true.)) error stop

end program example_math_is_close
