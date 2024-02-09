program example_eye1
  use stdlib_linalg, only: eye
  implicit none
  integer :: i(2, 2)
  real :: a(3, 3)
  real :: b(2, 3)  !! Matrix is non-square.
  complex :: c(2, 2)
  I = eye(2)              !! [1,0; 0,1]
  print *, "I = ", I
  print *, "sum(I) = ", sum(I)
  if (sum(I) /= 2) error stop
  if (I(0,0) /= 1 .and. I(1,1) /= 1) error stop
  A = eye(3)              !! [1.0,0.0,0.0; 0.0,1.0,0.0; 0.0,0.0,1.0]
  print *, "A = ", A
  print *, "sum(A) = ", sum(A)
  if (abs(sum(A) - 3.0) > 1e-8) error stop
  A = eye(3, 3)            !! [1.0,0.0,0.0; 0.0,1.0,0.0; 0.0,0.0,1.0]
  print *, "A = ", A
  print *, "sum(A) = ", sum(A)
  if (abs(sum(A) - 3.0) > 1e-8) error stop
  B = eye(2, 3)            !! [1.0,0.0,0.0; 0.0,1.0,0.0]
  print *, "B = ", B
  print *, "sum(B) = ", sum(B)
  if (abs(sum(B) - 2.0) > 1e-8) error stop
  C = eye(2, 2)            !! [(1.0,0.0),(0.0,0.0); (0.0,0.0),(1.0,0.0)]
  print *, "C = ", C
  print *, "abs(sum(C)) = ", abs(sum(C))
  if (abs(abs(sum(C)) - 2.0) > 1e-8) error stop
  C = eye(2, 2)
  C = C * cmplx(1.0, 1.0)!! [(1.0,1.0),(0.0,0.0); (0.0,0.0),(1.0,1.0)]
  print *, "C = ", C
  print *, "abs(sum(C)) = ", abs(sum(C))
  if (abs(abs(sum(C)) - 2.828427) > 1e-8) error stop
end program example_eye1
