program example_diff

  use stdlib_math, only: diff
  implicit none

  integer :: i(7) = [1, 1, 2, 3, 5, 8, 13]
  real    :: x(6) = [0, 5, 15, 30, 50, 75]
  integer :: A(3, 3) = reshape([1, 7, 17, 3, 11, 19, 5, 13, 23], [3, 3])
  integer :: Y(3, 2)

  print *, diff(i)        ! [0, 1, 1, 2, 3, 5]
  print *, sum(diff(i))
  if (sum(diff(i)) /= 12) error stop

  print *, diff(x, 2)     ! [5.0, 5.0, 5.0, 5.0]
  print *, sum(diff(x, 2))
  if (abs(sum(diff(x, 2)) - 20.0) > 1e-8) error stop

  Y = diff(A, n=1, dim=2)

  print *, Y(1, :)        ! [2, 2]
  print *, sum(Y(1, :))
  if (sum(Y(1, :)) /= 4) error stop
  if (any(Y(1, :) /= 2)) error stop

  print *, Y(2, :)        ! [4, 2]
  print *, sum(Y(2, :))
  if (sum(Y(2, :)) /= 6) error stop
  if (Y(2,1) /= 4) error stop
  if (Y(2,2) /= 2) error stop

  print *, Y(3, :)        ! [2, 4]
  print *, sum(Y(3, :))
  if (sum(Y(3, :)) /= 6) error stop
  if (Y(3,1) /= 2) error stop
  if (Y(3,2) /= 4) error stop

  print *, diff(i, prepend=[0]) ! [1, 0, 1, 1, 2, 3, 5]
  print *, sum(diff(i, prepend=[0]))
  if (sum(diff(i, prepend=[0])) /= 13) error stop

  print *, diff(i, append=[21]) ! [0, 1, 1, 2, 3, 5, 8]
  print *, sum(diff(i, append=[21]))
  if (sum(diff(i, append=[21])) /= 20) error stop

end program example_diff
