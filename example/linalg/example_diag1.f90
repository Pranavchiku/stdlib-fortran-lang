program example_diag1
  use stdlib_linalg, only: diag
  implicit none
  real, allocatable :: A(:, :)
  integer :: i
  allocate(A(10, 10))
  A = diag([(1, i=1, 10)]) ! creates a 10 by 10 identity matrix
end program example_diag1
