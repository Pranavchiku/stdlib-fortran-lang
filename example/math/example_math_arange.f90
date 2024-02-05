program example_math_arange
  use stdlib_math, only: arange

  print *, arange(3)                 ! [1,2,3]
  if (sum(arange(3)) /= 6) error stop

  print *, arange(-1)                ! [1,0,-1]
  ! if (sum(arange(-1)) /= 0) error stop ! TODO: fix this

  print *, arange(0, 2)               ! [0,1,2]
  if (sum(arange(0, 2)) /= 3) error stop


  print *, arange(1, -1)              ! [1,0,-1]
  ! if (sum(arange(1, -1)) /= 0) error stop ! TODO: fix this

  print *, arange(0, 2, 2)           ! [0,2]
  if (sum(arange(0, 2, 2)) /= 2) error stop

  print *, arange(3.0)               ! [1.0,2.0,3.0]
  if (abs(sum(arange(3.0)) - 6.0) > 1.0e-8) error stop

  print *, arange(0.0, 5.0)           ! [0.0,1.0,2.0,3.0,4.0,5.0]
  if (abs(sum(arange(0.0, 5.0)) - 15.0) > 1.0e-8) error stop

  print *, arange(0.0, 6.0, 2.5)       ! [0.0,2.5,5.0]
  if (abs(sum(arange(0.0, 6.0, 2.5)) - 7.5) > 1.0e-8) error stop

  ! print *, cmplx(1.0, 1.0)*arange(3)       ! [(1.0,1.0),(2.0,2.0),[3.0,3.0]]

  print *, arange(0.0, 2.0, -2.0)      ! [0.0,2.0].     Not recommended: `step` argument is negative!
  if (abs(sum(arange(0.0, 2.0, -2.0)) - 2.0) > 1.0e-8) error stop

  print *, arange(0.0, 2.0, 0.0)       ! [0.0,1.0,2.0]. Not recommended: `step` argument is zero!
  if (abs(sum(arange(0.0, 2.0, 0.0)) - 3.0) > 1.0e-8) error stop

end program example_math_arange
