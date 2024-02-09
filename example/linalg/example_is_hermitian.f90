program example_is_hermitian
  use stdlib_linalg, only: is_hermitian
  implicit none
  complex :: A(2, 2), B(2, 2)
  logical :: res
  complex, parameter :: c(4) = [cmplx(1., 0.), cmplx(3., -1.), cmplx(3., 1.), cmplx(4., 0.)]
  complex, parameter :: d(4) = [cmplx(1., 0.), cmplx(3., 1.), cmplx(3., 1.), cmplx(4., 0.)]
  A = reshape(c, shape(A))
  B = reshape(d, shape(B))
  res = is_hermitian(A) ! returns .true.
  print *, "res = ", res
  if (.not. res) error stop
  res = is_hermitian(B) ! returns .false.
  print *, "res = ", res
  ! if (res) error stop ! TODO: fix this
end program example_is_hermitian
