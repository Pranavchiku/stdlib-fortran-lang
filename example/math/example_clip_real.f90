program example_clip_real
  use stdlib_math, only: clip
  use stdlib_kinds, only: sp
  implicit none
  real(sp) :: x
  real(sp) :: xmin
  real(sp) :: xmax
  real(sp) :: clipped_value

  xmin = -5.769_sp
  xmax = 3.025_sp
  x = 3.025_sp

  clipped_value = clip(x, xmin, xmax)
  print *, "clipped_value <- ", clipped_value
  if (abs(clipped_value - 3.025_sp) > 1e-8) error stop
! clipped_value <- 3.02500010
end program example_clip_real
