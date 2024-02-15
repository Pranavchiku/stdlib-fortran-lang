program example_trace
  use stdlib_linalg, only: trace
  implicit none
  real :: A(3, 3)
  real, parameter :: x(9) = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  A = reshape(x, [3, 3])
  print *, trace(A) ! 1 + 5 + 9
  if (trace(A) /= 15) error stop
end program example_trace
