program example_math_arg
  use stdlib_math, only: arg
  real :: res(4)
  integer :: i
  complex :: z(4)
  print *, arg((0.0, 0.0))                  ! 0.0
  if (abs(arg((0.0, 0.0)) - 0.0) > 1.0e-8) error stop

  print *, arg((3.0, 4.0))                  ! 0.927
  if (abs(arg((3.0, 4.0)) - 0.9272952180016122) > 1.0e-8) error stop

  print *, arg(2.0*exp((0.0, 0.5)))         ! 0.5
  if (abs(arg(2.0*exp((0.0, 0.5))) - 0.5) > 1.0e-8) error stop

  print *, arg([(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)])  ! [π/2, 0.0, -π/2, π]

  z = [(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)]
  do i = 1, 4
    res(i) = arg(z(i))
    print *, res(i)
  end do

  if (abs(res(1) - 1.57079637) > 1.0e-8) error stop
  if (abs(res(2) - 0.0) > 1.0e-8) error stop
  if (abs(res(3) + 1.57079637) > 1.0e-8) error stop
  if (abs(res(4) - 3.1415927) > 1.0e-8) error stop

end program example_math_arg
