program example_diag5
  use stdlib_linalg, only: diag
  implicit none
  integer, parameter :: n = 3
  real :: A(n, n)
  real, allocatable :: v(:)
  real, parameter :: x(9) = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  A = reshape(x, [n, n])
  v = diag(A, -1) ! v is [2,6]
  v = diag(A, 1)  ! v is [4,8]
end program example_diag5
