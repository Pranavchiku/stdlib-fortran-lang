module stdlib_sorting_radix_sort

    use stdlib_kinds, only: &
        int8,               &
        int16,              &
        int32,              &
        int64,              &
        sp,                 &
        dp 
    
    use stdlib_optval, only: optval

    implicit none

    integer, parameter :: radix_bits = 8
    integer, parameter :: radix_mask = 255
    integer(kind=int16), parameter :: radix_bits_i16 = 8_int16
    integer(kind=int16), parameter :: radix_mask_i16 = 255_int16
    integer(kind=int32), parameter :: radix_bits_i32 = 8_int32
    integer(kind=int32), parameter :: radix_mask_i32 = 255_int32
    integer(kind=int64), parameter :: radix_bits_i64 = 8_int64
    integer(kind=int64), parameter :: radix_mask_i64 = 255_int64

    private

    integer, parameter, public :: int_size = int64 !! Integer kind for indexing

    public radix_sort

    interface radix_sort
        module procedure int8_radix_sort
        module procedure int16_radix_sort
        module procedure int32_radix_sort
        module procedure int64_radix_sort
        module procedure sp_radix_sort
        module procedure dp_radix_sort
    end interface radix_sort

contains
    pure subroutine int8_radix_sort(array, reverse)
        integer(kind=int8), dimension(:), intent(inout) :: array
        logical, intent(in), optional :: reverse
        integer(kind=int8) :: item
        integer(kind=int_size) :: i, N
        
        integer :: bin_idx
        integer(kind=int_size), dimension(-128:127) :: counts

        N = size(array, kind=int_size)

        counts(:) = 0
        do i = 1, N
            bin_idx = array(i)
            counts(bin_idx) = counts(bin_idx) + 1
        end do
        i = 1
        do bin_idx = -128, 127
            do while (counts(bin_idx) > 0)
                array(i) = int(bin_idx, kind=int8)
                i = i + 1
                counts(bin_idx) = counts(bin_idx) - 1
            end do
        end do

        if (optval(reverse, .false.)) then
            do i = 1, N/2
                item = array(i)
                array(i) = array(N - i + 1)
                array(N - i + 1) = item
            end do
        end if
    end subroutine int8_radix_sort

    pure subroutine int16_radix_sort(array, work, reverse)
        integer(kind=int16), dimension(:), intent(inout) :: array
        integer(kind=int16), dimension(:), intent(inout), target, optional :: work
        logical, intent(in), optional :: reverse
        integer(kind=int_size) :: i, N, start, middle, end
        integer(kind=int16), dimension(:), pointer :: buffer
        integer(kind=int16) :: item
        logical :: use_internal_buffer
        
        integer :: b, b0, b1
        integer(kind=int_size), dimension(0:radix_mask) :: c0, c1

        N = size(array, kind=int_size)
        if (present(work)) then
            if (size(work, kind=int_size) < N) then
                error stop "int16_radix_sort: work array is too small."
            end if
            use_internal_buffer = .false.
            buffer => work
        else
            use_internal_buffer = .true.
            allocate (buffer(N))
        end if

        c0(:) = 0
        c1(:) = 0
        do i = 1, N
            b0 = iand(array(i), radix_mask_i16)
            b1 = ishft(array(i), -radix_bits_i16)
            c0(b0) = c0(b0) + 1
            c1(b1) = c1(b1) + 1
        end do
        do b = 1, radix_mask
            c0(b) = c0(b) + c0(b - 1)
            c1(b) = c1(b) + c1(b - 1)
        end do
        do i = N, 1, -1
            b0 = iand(array(i), radix_mask_i16)
            buffer(c0(b0)) = array(i)
            c0(b0) = c0(b0) - 1
        end do
        do i = N, 1, -1
            b1 = ishft(buffer(i), -radix_bits_i16)
            array(c1(b1)) = buffer(i)
            c1(b1) = c1(b1) - 1
        end do

! After calling `radix_sort_u<width>_helper. The array is sorted as unsigned integers.
! In the view of signed array, the negative numbers are sorted but in the tail of the array.
! Use binary search to find the boundary, and move them to correct position.
        if (array(1) >= 0 .and. array(N) < 0) then
            start = 1
            end = N
            middle = (1 + N)/2
            do while (.true.)
                if (array(middle) >= 0) then
                    start = middle + 1
                else
                    end = middle
                end if
                middle = (start + end)/2
                if (start == end) exit
            end do
            buffer(1:(N - middle + 1)) = array(middle:N)
            buffer(N - middle + 2:N) = array(1:middle - 1)
            array(:) = buffer(:)
        end if
        if (optval(reverse, .false.)) then
            do i = 1, N/2
                item = array(i)
                array(i) = array(N - i + 1)
                array(N - i + 1) = item
            end do
        end if
        if (use_internal_buffer) then
            deallocate (buffer)
        end if
    end subroutine int16_radix_sort

    pure subroutine int32_radix_sort(array, work, reverse)
        integer(kind=int32), dimension(:), intent(inout) :: array
        integer(kind=int32), dimension(:), intent(inout), target, optional :: work
        logical, intent(in), optional :: reverse
        integer(kind=int_size) :: i, N, start, middle, end
        integer(kind=int32), dimension(:), pointer :: buffer
        integer(kind=int32) :: item
        logical :: use_internal_buffer

        integer :: b, b0, b1, b2, b3
        integer(kind=int_size), dimension(0:radix_mask) :: c0, c1, c2, c3

        N = size(array, kind=int_size)
        if (present(work)) then
            if (size(work, kind=int_size) < N) then
                error stop "int32_radix_sort: work array is too small."
            end if
            use_internal_buffer = .false.
            buffer => work
        else
            use_internal_buffer = .true.
            allocate (buffer(N))
        end if

        c0(:) = 0
        c1(:) = 0
        c2(:) = 0
        c3(:) = 0
        do i = 1, N
            b0 = iand(array(i), radix_mask_i32)
            b1 = iand(ishft(array(i), -radix_bits_i32), radix_mask_i32)
            b2 = iand(ishft(array(i), -2*radix_bits_i32), radix_mask_i32)
            b3 = ishft(array(i), -3*radix_bits_i32)
            c0(b0) = c0(b0) + 1
            c1(b1) = c1(b1) + 1
            c2(b2) = c2(b2) + 1
            c3(b3) = c3(b3) + 1
        end do
        do b = 1, radix_mask
            c0(b) = c0(b) + c0(b - 1)
            c1(b) = c1(b) + c1(b - 1)
            c2(b) = c2(b) + c2(b - 1)
            c3(b) = c3(b) + c3(b - 1)
        end do
        do i = N, 1, -1
            b0 = iand(array(i), radix_mask_i32)
            buffer(c0(b0)) = array(i)
            c0(b0) = c0(b0) - 1
        end do
        do i = N, 1, -1
            b1 = iand(ishft(buffer(i), -radix_bits_i32), radix_mask_i32)
            array(c1(b1)) = buffer(i)
            c1(b1) = c1(b1) - 1
        end do
        do i = N, 1, -1
            b2 = iand(ishft(array(i), -2*radix_bits_i32), radix_mask_i32)
            buffer(c2(b2)) = array(i)
            c2(b2) = c2(b2) - 1
        end do
        do i = N, 1, -1
            b3 = ishft(buffer(i), -3*radix_bits_i32)
            array(c3(b3)) = buffer(i)
            c3(b3) = c3(b3) - 1
        end do

        if (array(1) >= 0 .and. array(N) < 0) then
            start = 1
            end = N
            middle = (1 + N)/2
            do while (.true.)
                if (array(middle) >= 0) then
                    start = middle + 1
                else
                    end = middle
                end if
                middle = (start + end)/2
                if (start == end) exit
            end do
            buffer(1:(N - middle + 1)) = array(middle:N)
            buffer(N - middle + 2:N) = array(1:middle - 1)
            array(:) = buffer(:)
        end if
        if (optval(reverse, .false.)) then
            do i = 1, N/2
                item = array(i)
                array(i) = array(N - i + 1)
                array(N - i + 1) = item
            end do
        end if
        if (use_internal_buffer) then
            deallocate (buffer)
        end if
    end subroutine int32_radix_sort

    pure subroutine sp_radix_sort(array, work, reverse)
        use iso_c_binding
        real(kind=sp), dimension(:), intent(inout), target :: array
        real(kind=sp), dimension(:), intent(inout), target, optional :: work
        logical, intent(in), optional :: reverse
        integer(kind=int_size) :: i, N, pos, rev_pos
        integer(kind=int32), dimension(:), pointer :: arri32
        integer(kind=int32), dimension(:), pointer :: buffer
        real(kind=sp) :: item
        logical :: use_internal_buffer

        integer :: b, b0, b1, b2, b3
        integer(kind=int_size), dimension(0:radix_mask) :: c0, c1, c2, c3

        N = size(array, kind=int_size)
        if (present(work)) then
            if (size(work, kind=int_size) < N) then
                error stop "sp_radix_sort: work array is too small."
            end if
            use_internal_buffer = .false.
            call c_f_pointer(c_loc(work), buffer, [N])
        else
            use_internal_buffer = .true.
            allocate (buffer(N))
        end if
        call c_f_pointer(c_loc(array), arri32, [N])

        c0(:) = 0
        c1(:) = 0
        c2(:) = 0
        c3(:) = 0
        do i = 1, N
            b0 = iand(arri32(i), radix_mask_i32)
            b1 = iand(ishft(arri32(i), -radix_bits_i32), radix_mask_i32)
            b2 = iand(ishft(arri32(i), -2*radix_bits_i32), radix_mask_i32)
            b3 = ishft(arri32(i), -3*radix_bits_i32)
            c0(b0) = c0(b0) + 1
            c1(b1) = c1(b1) + 1
            c2(b2) = c2(b2) + 1
            c3(b3) = c3(b3) + 1
        end do
        do b = 1, radix_mask
            c0(b) = c0(b) + c0(b - 1)
            c1(b) = c1(b) + c1(b - 1)
            c2(b) = c2(b) + c2(b - 1)
            c3(b) = c3(b) + c3(b - 1)
        end do
        do i = N, 1, -1
            b0 = iand(arri32(i), radix_mask_i32)
            buffer(c0(b0)) = arri32(i)
            c0(b0) = c0(b0) - 1
        end do
        do i = N, 1, -1
            b1 = iand(ishft(buffer(i), -radix_bits_i32), radix_mask_i32)
            arri32(c1(b1)) = buffer(i)
            c1(b1) = c1(b1) - 1
        end do
        do i = N, 1, -1
            b2 = iand(ishft(arri32(i), -2*radix_bits_i32), radix_mask_i32)
            buffer(c2(b2)) = arri32(i)
            c2(b2) = c2(b2) - 1
        end do
        do i = N, 1, -1
            b3 = ishft(buffer(i), -3*radix_bits_i32)
            arri32(c3(b3)) = buffer(i)
            c3(b3) = c3(b3) - 1
        end do

! After calling `radix_sort_u<width>_helper. The array is sorted as unsigned integers.
! The positive real number is sorted, guaranteed by IEEE-754 standard.
! But the negative real number is sorted in a reversed order, and also in the tail of array.
! Remark that -0.0 is the minimum nagative integer, so using radix sort, -0.0 is naturally lesser than 0.0.
! In IEEE-754 standard, the bit representation of `Inf` is greater than all positive real numbers,
! and the `-Inf` is lesser than all negative real numbers. So the order of `Inf, -Inf` is ok.
! The bit representation of `NaN` may be positive or negative integers in different machine,
! thus if the array contains `NaN`, the result is undefined.
        if (arri32(1) >= 0 .and. arri32(N) < 0) then
            pos = 1
            rev_pos = N
            do while (arri32(rev_pos) < 0)
                buffer(pos) = arri32(rev_pos)
                pos = pos + 1
                rev_pos = rev_pos - 1
            end do
            buffer(pos:N) = arri32(1:rev_pos)
            arri32(:) = buffer(:)
        end if
        if (optval(reverse, .false.)) then
            do i = 1, N/2
                item = array(i)
                array(i) = array(N - i + 1)
                array(N - i + 1) = item
            end do
        end if
        if (use_internal_buffer) then
            deallocate (buffer)
        end if
    end subroutine sp_radix_sort

    pure subroutine int64_radix_sort(array, work, reverse)
        integer(kind=int64), dimension(:), intent(inout) :: array
        integer(kind=int64), dimension(:), intent(inout), target, optional :: work
        logical, intent(in), optional :: reverse
        integer(kind=int_size) :: i, N, start, middle, end
        integer(kind=int64), dimension(:), pointer :: buffer
        integer(kind=int64) :: item
        logical :: use_internal_buffer

        integer(kind=int64) :: b, b0, b1, b2, b3, b4, b5, b6, b7
        integer(kind=int_size), dimension(0:radix_mask) :: c0, c1, c2, c3, c4, c5, c6, c7

        N = size(array, kind=int_size)
        if (present(work)) then
            if (size(work, kind=int_size) < N) then
                error stop "int64_radix_sort: work array is too small."
            end if
            use_internal_buffer = .false.
            buffer => work
        else
            use_internal_buffer = .true.
            allocate (buffer(N))
        end if
        
        c0(:) = 0
        c1(:) = 0
        c2(:) = 0
        c3(:) = 0
        c4(:) = 0
        c5(:) = 0
        c6(:) = 0
        c7(:) = 0
        do i = 1, N
            b0 = iand(array(i), radix_mask_i64)
            b1 = iand(ishft(array(i), -radix_bits_i64), radix_mask_i64)
            b2 = iand(ishft(array(i), -2*radix_bits_i64), radix_mask_i64)
            b3 = iand(ishft(array(i), -3*radix_bits_i64), radix_mask_i64)
            b4 = iand(ishft(array(i), -4*radix_bits_i64), radix_mask_i64)
            b5 = iand(ishft(array(i), -5*radix_bits_i64), radix_mask_i64)
            b6 = iand(ishft(array(i), -6*radix_bits_i64), radix_mask_i64)
            b7 = ishft(array(i), -7*radix_bits_i64)
            c0(b0) = c0(b0) + 1
            c1(b1) = c1(b1) + 1
            c2(b2) = c2(b2) + 1
            c3(b3) = c3(b3) + 1
            c4(b4) = c4(b4) + 1
            c5(b5) = c5(b5) + 1
            c6(b6) = c6(b6) + 1
            c7(b7) = c7(b7) + 1
        end do
        do b = 1, radix_mask
            c0(b) = c0(b) + c0(b - 1)
            c1(b) = c1(b) + c1(b - 1)
            c2(b) = c2(b) + c2(b - 1)
            c3(b) = c3(b) + c3(b - 1)
            c4(b) = c4(b) + c4(b - 1)
            c5(b) = c5(b) + c5(b - 1)
            c6(b) = c6(b) + c6(b - 1)
            c7(b) = c7(b) + c7(b - 1)
        end do
        do i = N, 1, -1
            b0 = iand(array(i), radix_mask_i64)
            buffer(c0(b0)) = array(i)
            c0(b0) = c0(b0) - 1
        end do
        do i = N, 1, -1
            b1 = iand(ishft(buffer(i), -radix_bits_i64), radix_mask_i64)
            array(c1(b1)) = buffer(i)
            c1(b1) = c1(b1) - 1
        end do
        do i = N, 1, -1
            b2 = iand(ishft(array(i), -2*radix_bits_i64), radix_mask_i64)
            buffer(c2(b2)) = array(i)
            c2(b2) = c2(b2) - 1
        end do
        do i = N, 1, -1
            b3 = iand(ishft(buffer(i), -3*radix_bits_i64), radix_mask_i64)
            array(c3(b3)) = buffer(i)
            c3(b3) = c3(b3) - 1
        end do
        do i = N, 1, -1
            b4 = iand(ishft(array(i), -4*radix_bits_i64), radix_mask_i64)
            buffer(c4(b4)) = array(i)
            c4(b4) = c4(b4) - 1
        end do
        do i = N, 1, -1
            b5 = iand(ishft(buffer(i), -5*radix_bits_i64), radix_mask_i64)
            array(c5(b5)) = buffer(i)
            c5(b5) = c5(b5) - 1
        end do
        do i = N, 1, -1
            b6 = iand(ishft(array(i), -6*radix_bits_i64), radix_mask_i64)
            buffer(c6(b6)) = array(i)
            c6(b6) = c6(b6) - 1
        end do
        do i = N, 1, -1
            b7 = ishft(buffer(i), -7*radix_bits_i64)
            array(c7(b7)) = buffer(i)
            c7(b7) = c7(b7) - 1
        end do
        
        if (array(1) >= 0 .and. array(N) < 0) then
            start = 1
            end = N
            middle = (1 + N)/2
            do while (.true.)
                if (array(middle) >= 0) then
                    start = middle + 1
                else
                    end = middle
                end if
                middle = (start + end)/2
                if (start == end) exit
            end do
            buffer(1:(N - middle + 1)) = array(middle:N)
            buffer(N - middle + 2:N) = array(1:middle - 1)
            array(:) = buffer(:)
        end if
        if (optval(reverse, .false.)) then
            do i = 1, N/2
                item = array(i)
                array(i) = array(N - i + 1)
                array(N - i + 1) = item
            end do
        end if
        if (use_internal_buffer) then
            deallocate (buffer)
        end if
    end subroutine int64_radix_sort

    subroutine dp_radix_sort(array, work, reverse)
        use iso_c_binding
        real(kind=dp), dimension(:), intent(inout), target :: array
        real(kind=dp), dimension(:), intent(inout), target, optional :: work
        logical, intent(in), optional :: reverse
        integer(kind=int_size) :: i, N, pos, rev_pos
        integer(kind=int64), dimension(:), pointer :: arri64
        integer(kind=int64), dimension(:), pointer :: buffer
        real(kind=dp) :: item
        logical :: use_internal_buffer

        integer(kind=int64) :: b, b0, b1, b2, b3, b4, b5, b6, b7
        integer(kind=int_size), dimension(0:radix_mask) :: c0, c1, c2, c3, c4, c5, c6, c7

        N = size(array, kind=int_size)
        if (present(work)) then
            if (size(work, kind=int_size) < N) then
                error stop "sp_radix_sort: work array is too small."
            end if
            use_internal_buffer = .false.
            call c_f_pointer(c_loc(work), buffer, [N])
        else
            use_internal_buffer = .true.
            allocate (buffer(N))
        end if
        call c_f_pointer(c_loc(array), arri64, [N])

        c0(:) = 0
        c1(:) = 0
        c2(:) = 0
        c3(:) = 0
        c4(:) = 0
        c5(:) = 0
        c6(:) = 0
        c7(:) = 0
        do i = 1, N
            b0 = iand(arri64(i), radix_mask_i64)
            b1 = iand(ishft(arri64(i), -radix_bits_i64), radix_mask_i64)
            b2 = iand(ishft(arri64(i), -2*radix_bits_i64), radix_mask_i64)
            b3 = iand(ishft(arri64(i), -3*radix_bits_i64), radix_mask_i64)
            b4 = iand(ishft(arri64(i), -4*radix_bits_i64), radix_mask_i64)
            b5 = iand(ishft(arri64(i), -5*radix_bits_i64), radix_mask_i64)
            b6 = iand(ishft(arri64(i), -6*radix_bits_i64), radix_mask_i64)
            b7 = ishft(arri64(i), -7*radix_bits_i64)
            c0(b0) = c0(b0) + 1
            c1(b1) = c1(b1) + 1
            c2(b2) = c2(b2) + 1
            c3(b3) = c3(b3) + 1
            c4(b4) = c4(b4) + 1
            c5(b5) = c5(b5) + 1
            c6(b6) = c6(b6) + 1
            c7(b7) = c7(b7) + 1
        end do
        do b = 1, radix_mask
            c0(b) = c0(b) + c0(b - 1)
            c1(b) = c1(b) + c1(b - 1)
            c2(b) = c2(b) + c2(b - 1)
            c3(b) = c3(b) + c3(b - 1)
            c4(b) = c4(b) + c4(b - 1)
            c5(b) = c5(b) + c5(b - 1)
            c6(b) = c6(b) + c6(b - 1)
            c7(b) = c7(b) + c7(b - 1)
        end do
        do i = N, 1, -1
            b0 = iand(arri64(i), radix_mask_i64)
            buffer(c0(b0)) = arri64(i)
            c0(b0) = c0(b0) - 1
        end do
        do i = N, 1, -1
            b1 = iand(ishft(buffer(i), -radix_bits_i64), radix_mask_i64)
            arri64(c1(b1)) = buffer(i)
            c1(b1) = c1(b1) - 1
        end do
        do i = N, 1, -1
            b2 = iand(ishft(arri64(i), -2*radix_bits_i64), radix_mask_i64)
            buffer(c2(b2)) = arri64(i)
            c2(b2) = c2(b2) - 1
        end do
        do i = N, 1, -1
            b3 = iand(ishft(buffer(i), -3*radix_bits_i64), radix_mask_i64)
            arri64(c3(b3)) = buffer(i)
            c3(b3) = c3(b3) - 1
        end do
        do i = N, 1, -1
            b4 = iand(ishft(arri64(i), -4*radix_bits_i64), radix_mask_i64)
            buffer(c4(b4)) = arri64(i)
            c4(b4) = c4(b4) - 1
        end do
        do i = N, 1, -1
            b5 = iand(ishft(buffer(i), -5*radix_bits_i64), radix_mask_i64)
            arri64(c5(b5)) = buffer(i)
            c5(b5) = c5(b5) - 1
        end do
        do i = N, 1, -1
            b6 = iand(ishft(arri64(i), -6*radix_bits_i64), radix_mask_i64)
            buffer(c6(b6)) = arri64(i)
            c6(b6) = c6(b6) - 1
        end do
        do i = N, 1, -1
            b7 = ishft(buffer(i), -7*radix_bits_i64)
            arri64(c7(b7)) = buffer(i)
            c7(b7) = c7(b7) - 1
        end do

        if (arri64(1) >= 0 .and. arri64(N) < 0) then
            pos = 1
            rev_pos = N
            do while (arri64(rev_pos) < 0)
                buffer(pos) = arri64(rev_pos)
                pos = pos + 1
                rev_pos = rev_pos - 1
            end do
            buffer(pos:N) = arri64(1:rev_pos)
            arri64(:) = buffer(:)
        end if
        if (optval(reverse, .false.)) then
            do i = 1, N/2
                item = array(i)
                array(i) = array(N - i + 1)
                array(N - i + 1) = item
            end do
        end if
        if (use_internal_buffer) then
            deallocate (buffer)
        end if
    end subroutine dp_radix_sort
end module stdlib_sorting_radix_sort
