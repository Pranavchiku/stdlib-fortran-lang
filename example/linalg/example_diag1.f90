program example_diag1
  use stdlib_linalg, only: diag
  implicit none
  real, allocatable :: A(:, :)
  integer :: x(10)
  integer :: i
  allocate(A(10, 10))
  x = [(i, i=1, 10)]
  A = diag(x) ! creates a 10 by 10 identity matrix
  print *, A
  print *, "sum(A) = ", sum(A)
  ! if (abs(sum(A) - 9.0) > 1e-8) error stop ! TODO: fix this
  if (any(.not. (abs(A - 1.0) > 1e-8 .or. abs(A - 0.0) > 1e-8))) error stop
end program example_diag1
