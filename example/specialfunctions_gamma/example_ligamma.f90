program example_ligamma
  use stdlib_specialfunctions_gamma, only: lig => lower_incomplete_gamma
  implicit none
  integer :: p
  real :: p1

  p = 3
  p1 = 2.3
  print *, lig(p, -5.0)
  if (abs(lig(p, -5.0) + 2521.02417) > 1e-8) error stop

! -2521.02417

  print *, lig(p1, 5.0)
  if (abs(lig(p1, 5.0) - 1.09715652) > 1e-8) error stop

! 1.09715652
end program example_ligamma
