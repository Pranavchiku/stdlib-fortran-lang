program example_fwrite
  use stdlib_string_type
  implicit none
  type(string_type) :: string
  integer :: io, iostat
  character(len=100) :: iomsg
  string = "Important saved value"

  open (newunit=io, form="formatted", status="scratch")
  write (io, *, iostat=iostat, iomsg=iomsg) string
  write (io, *, iostat=iostat, iomsg=iomsg)

  rewind (io)

  read (io, *, iostat=iostat, iomsg=iomsg) string
  close (io)
end program example_fwrite
