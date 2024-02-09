program example_math_argpi
  use stdlib_math, only: argpi
  integer :: i
  real :: res(4)
  complex :: z(4)
  print *, argpi((0.0, 0.0))                  ! 0.0
  if (abs(argpi((0.0, 0.0)) - 0.0) > 1.0e-8) error stop

  print *, argpi((3.0, 4.0))                  ! 0.295
  if (abs(argpi((3.0, 4.0)) - 0.295167238) > 1.0e-8) error stop
  print *, argpi(2.0*exp((0.0, 0.5)))         ! 0.159
  if (abs(argpi(2.0*exp((0.0, 0.5))) - 0.159154937) > 1.0e-8) error stop
  print *, argpi([(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)])  ! [0.5, 0.0, -0.5, 1.0]
  z = [(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)]
  do i = 1, 4
    res(i) = argpi(z(i))
  end do
  if (abs(res(1) - 0.5) > 1.0e-8) error stop
  if (abs(res(2) - 0.0) > 1.0e-8) error stop
  if (abs(res(3) + 0.5) > 1.0e-8) error stop
  if (abs(res(4) - 1.0) > 1.0e-8) error stop
end program example_math_argpi
