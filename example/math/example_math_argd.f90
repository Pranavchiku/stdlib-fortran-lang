program example_math_argd
  use stdlib_math, only: argd
  integer :: i
  real :: res(4)
  complex :: z(4)
  print *, argd((0.0, 0.0))                  ! 0.0°
  if (abs(argd((0.0, 0.0)) - 0.0) > 1e-8) error stop

  print *, argd((3.0, 4.0))                  ! 53.1°
  if (abs(argd((3.0, 4.0)) - 53.1301003) > 1e-8) error stop

  print *, argd(2.0*exp((0.0, 0.5)))         ! 28.64°
  if (abs(argd(2.0*exp((0.0, 0.5))) - 28.6478882) > 1e-8) error stop

  print *, argd([(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)])  ! [90°, 0°, -90°, 180°]
  z = [(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)]
  do i = 1, 4
    res(i) = argd(z(i))
  end do

  if (abs(res(1) - 90.0) > 1e-8) error stop
  if (abs(res(2) - 0.0) > 1e-8) error stop
  if (abs(res(3) + 90.0) > 1e-8) error stop
  if (abs(res(4) - 179.999985) > 1e-8) error stop

end program example_math_argd
