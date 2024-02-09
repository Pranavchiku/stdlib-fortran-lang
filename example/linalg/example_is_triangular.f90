program example_is_triangular
  use stdlib_linalg, only: is_triangular
  implicit none
  real :: A(3, 3), B(3, 3)
  logical :: res
  A = reshape([1., 0., 0., 4., 5., 0., 7., 8., 9.], shape(A))
  B = reshape([1., 0., 3., 4., 5., 0., 7., 8., 9.], shape(B))
  res = is_triangular(A, 'u') ! returns .true.
  if (.not. res) error stop
  res = is_triangular(B, 'u') ! returns .false.
  if (res) error stop
end program example_is_triangular
