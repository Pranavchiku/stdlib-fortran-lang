program example_eye2
  use stdlib_linalg, only: eye, diag
  implicit none
  print *, all(eye(4) == diag([1, 1, 1, 1])) ! prints .true.
  if (.not. all(eye(4) == diag([1, 1, 1, 1]))) error stop
end program example_eye2
