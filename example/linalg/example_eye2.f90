program example_eye2
  use stdlib_linalg, only: eye, diag
  implicit none
  ! print *, all(eye(4) == diag([1, 1, 1, 1])) ! prints .true. ! TODO: Fix this
  real :: A(4, 4)
  real :: B(4, 4)
  A = eye(4)
  print *, A
  B = diag([1, 1, 1, 1])
  print *, B
  print *, all(A == B) ! prints .true.
end program example_eye2
