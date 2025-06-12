program example_uread
  use stdlib_string_type
  implicit none
  type(string_type) :: string
  integer :: io, iostat
  character(len=100) :: iomsg
  string = "Important saved value"

  open (newunit=io, form="unformatted", status="scratch")
  write (io, iostat=iostat, iomsg=iomsg) string

  rewind (io)

  read (io, iostat=iostat, iomsg=iomsg) string
  close (io)
end program example_uread
