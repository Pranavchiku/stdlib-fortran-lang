program example_diag2
  use stdlib_linalg, only: diag
  implicit none
  real, allocatable :: v(:)
  real, allocatable :: A(:, :)
  v = [1, 2, 3, 4, 5]
  A = diag(v) ! creates a 5 by 5 matrix with elements of v on the diagonal
  print *, "A = ", A
  print *, "sum(A) = ", sum(A)
  if (abs(sum(A) - 15.0) > 1e-8) error stop
  if (abs(A(1,1) - 1.0) > 1e-8) error stop
  if (abs(A(2,2) - 2.0) > 1e-8) error stop
  if (abs(A(3,3) - 3.0) > 1e-8) error stop
  if (abs(A(4,4) - 4.0) > 1e-8) error stop
  if (abs(A(5,5) - 5.0) > 1e-8) error stop
end program example_diag2
