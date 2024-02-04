
module stdlib_math
    use stdlib_kinds, only: int8, int16, int32, int64, sp, dp, xdp, qp
    use stdlib_optval, only: optval

    implicit none
    private
    public :: clip, gcd, linspace, logspace
    public :: EULERS_NUMBER_SP, EULERS_NUMBER_DP
    public :: DEFAULT_LINSPACE_LENGTH, DEFAULT_LOGSPACE_BASE, DEFAULT_LOGSPACE_LENGTH
    public :: arange, arg, argd, argpi, is_close, all_close, diff

    integer, parameter :: DEFAULT_LINSPACE_LENGTH = 100
    integer, parameter :: DEFAULT_LOGSPACE_LENGTH = 50
    integer, parameter :: DEFAULT_LOGSPACE_BASE = 10

    real(sp), parameter :: sqrt_eps_sp = sqrt(epsilon(1.0_sp))
    real(dp), parameter :: sqrt_eps_dp = sqrt(epsilon(1.0_dp))

    ! Useful constants for lnspace
    real(sp), parameter :: EULERS_NUMBER_SP = exp(1.0_sp)
    real(dp), parameter :: EULERS_NUMBER_DP = exp(1.0_dp)

    !> Useful constants `PI` for `argd/argpi`
    real(kind=sp), parameter :: PI_sp = acos(-1.0_sp)
    real(kind=dp), parameter :: PI_dp = acos(-1.0_dp)

    interface clip
        module procedure clip_int8
        module procedure clip_int16
        module procedure clip_int32
        module procedure clip_int64
        module procedure clip_sp
        module procedure clip_dp
    end interface clip

    !> Returns the greatest common divisor of two integers
    !> ([Specification](../page/specs/stdlib_math.html#gcd))
    !>
    !> Version: experimental
    interface gcd
        module procedure gcd_int8
        module procedure gcd_int16
        module procedure gcd_int32
        module procedure gcd_int64
    end interface gcd

    interface linspace
    !! Version: Experimental
    !!
    !! Create rank 1 array of linearly spaced elements
    !! If the number of elements is not specified, create an array with size 100. If n is a negative value,
    !! return an array with size 0. If n = 1, return an array whose only element is end
    !!([Specification](../page/specs/stdlib_math.html#linspace-create-a-linearly-spaced-rank-one-array))
      procedure :: linspace_default_1_rsp_rsp
      procedure :: linspace_default_1_rdp_rdp
      procedure :: linspace_default_1_csp_csp
      procedure :: linspace_default_1_cdp_cdp

      procedure :: linspace_n_1_rsp_rsp
      procedure :: linspace_n_1_rdp_rdp

      procedure :: linspace_n_1_csp_csp
      procedure :: linspace_n_1_cdp_cdp

    ! Add support for integer linspace
    !!
    !! When dealing with integers as the `start` and `end` parameters, the return type is always a `real(dp)`.
      procedure :: linspace_default_1_iint8_iint8
      procedure :: linspace_default_1_iint16_iint16
      procedure :: linspace_default_1_iint32_iint32
      procedure :: linspace_default_1_iint64_iint64

      procedure :: linspace_n_1_iint8_iint8
      procedure :: linspace_n_1_iint16_iint16
      procedure :: linspace_n_1_iint32_iint32
      procedure :: linspace_n_1_iint64_iint64

  end interface

  interface logspace
  !! Version: Experimental
  !!
  !! Create rank 1 array of logarithmically spaced elements from base**start to base**end.
  !! If the number of elements is not specified, create an array with size 50. If n is a negative value,
  !! return an array with size 0. If n = 1, return an array whose only element is base**end. If no base
  !! is specified, logspace will default to using a base of 10
  !!
  !!([Specification](../page/specs/stdlib_math.html#logspace-create-a-logarithmically-spaced-rank-one-array))
    procedure :: logspace_1_rsp_default
    procedure :: logspace_1_rdp_default
    procedure :: logspace_1_csp_default
    procedure :: logspace_1_cdp_default
    procedure :: logspace_1_iint32_default

    procedure :: logspace_1_rsp_n
    procedure :: logspace_1_rdp_n
    procedure :: logspace_1_csp_n
    procedure :: logspace_1_cdp_n
    procedure :: logspace_1_iint32_n

    ! Generate logarithmically spaced sequence from sp base to the powers
    ! of sp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    procedure :: logspace_1_rsp_n_rbase

    procedure :: logspace_1_rsp_n_cbase

    procedure :: logspace_1_rsp_n_ibase
    ! Generate logarithmically spaced sequence from dp base to the powers
    ! of dp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    procedure :: logspace_1_rdp_n_rbase

    procedure :: logspace_1_rdp_n_cbase

    procedure :: logspace_1_rdp_n_ibase
    ! Generate logarithmically spaced sequence from sp base to the powers
    ! of sp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    procedure :: logspace_1_csp_n_rbase

    procedure :: logspace_1_csp_n_cbase

    procedure :: logspace_1_csp_n_ibase
    ! Generate logarithmically spaced sequence from dp base to the powers
    ! of dp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    procedure :: logspace_1_cdp_n_rbase

    procedure :: logspace_1_cdp_n_cbase

    procedure :: logspace_1_cdp_n_ibase
    ! Generate logarithmically spaced sequence from dp base to the powers
    ! of dp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    procedure :: logspace_1_iint32_n_rspbase

    procedure :: logspace_1_iint32_n_cspbase
    procedure :: logspace_1_iint32_n_rdpbase

    procedure :: logspace_1_iint32_n_cdpbase

    procedure :: logspace_1_iint32_n_ibase


  end interface

    !> Version: experimental
    !>
    !> `arange` creates a one-dimensional `array` of the `integer/real` type 
    !>  with fixed-spaced values of given spacing, within a given interval.
    !> ([Specification](../page/specs/stdlib_math.html#arange-function))
    interface arange
        procedure :: arange_r_sp
        procedure :: arange_r_dp
        procedure :: arange_i_int8
        procedure :: arange_i_int16
        procedure :: arange_i_int32
        procedure :: arange_i_int64
    end interface arange

    !> Version: experimental
    !>
    !> `arg` computes the phase angle in the interval (-π,π].
    !> ([Specification](../page/specs/stdlib_math.html#arg-function))
    interface arg
        procedure :: arg_sp
        procedure :: arg_dp
    end interface arg

    !> Version: experimental
    !>
    !> `argd` computes the phase angle of degree version in the interval (-180.0,180.0].
    !> ([Specification](../page/specs/stdlib_math.html#argd-function))
    interface argd
        procedure :: argd_sp
        procedure :: argd_dp
    end interface argd

    !> Version: experimental
    !>
    !> `argpi` computes the phase angle of circular version in the interval (-1.0,1.0].
    !> ([Specification](../page/specs/stdlib_math.html#argpi-function))
    interface argpi
        procedure :: argpi_sp
        procedure :: argpi_dp
    end interface argpi
    
    !> Returns a boolean scalar/array where two scalar/arrays are element-wise equal within a tolerance.
    !> ([Specification](../page/specs/stdlib_math.html#is_close-function))
    interface is_close
        procedure :: is_close_rsp
        procedure :: is_close_rdp
        procedure :: is_close_csp
        procedure :: is_close_cdp
    end interface is_close

    !> Version: experimental
    !>
    !> Returns a boolean scalar where two arrays are element-wise equal within a tolerance.
    !> ([Specification](../page/specs/stdlib_math.html#all_close-function))
    interface all_close
        procedure :: all_close_1_rsp
        procedure :: all_close_2_rsp
        procedure :: all_close_3_rsp
        procedure :: all_close_4_rsp
        procedure :: all_close_1_rdp
        procedure :: all_close_2_rdp
        procedure :: all_close_3_rdp
        procedure :: all_close_4_rdp
        procedure :: all_close_1_csp
        procedure :: all_close_2_csp
        procedure :: all_close_3_csp
        procedure :: all_close_4_csp
        procedure :: all_close_1_cdp
        procedure :: all_close_2_cdp
        procedure :: all_close_3_cdp
        procedure :: all_close_4_cdp
    end interface all_close
    
    !> Version: experimental
    !>
    !> Computes differences between adjacent elements of an array.
    !> ([Specification](../page/specs/stdlib_math.html#diff-function))
    interface diff
        procedure :: diff_1_sp
        procedure :: diff_2_sp
        procedure :: diff_1_dp
        procedure :: diff_2_dp
        procedure :: diff_1_int8
        procedure :: diff_2_int8
        procedure :: diff_1_int16
        procedure :: diff_2_int16
        procedure :: diff_1_int32
        procedure :: diff_2_int32
        procedure :: diff_1_int64
        procedure :: diff_2_int64
    end interface diff

contains

    elemental function clip_int8(x, xmin, xmax) result(res)
        integer(int8), intent(in) :: x
        integer(int8), intent(in) :: xmin
        integer(int8), intent(in) :: xmax
        integer(int8) :: res

        res = max(min(x, xmax), xmin)
    end function clip_int8

    elemental function clip_int16(x, xmin, xmax) result(res)
        integer(int16), intent(in) :: x
        integer(int16), intent(in) :: xmin
        integer(int16), intent(in) :: xmax
        integer(int16) :: res

        res = max(min(x, xmax), xmin)
    end function clip_int16

    elemental function clip_int32(x, xmin, xmax) result(res)
        integer(int32), intent(in) :: x
        integer(int32), intent(in) :: xmin
        integer(int32), intent(in) :: xmax
        integer(int32) :: res

        res = max(min(x, xmax), xmin)
    end function clip_int32

    elemental function clip_int64(x, xmin, xmax) result(res)
        integer(int64), intent(in) :: x
        integer(int64), intent(in) :: xmin
        integer(int64), intent(in) :: xmax
        integer(int64) :: res

        res = max(min(x, xmax), xmin)
    end function clip_int64

    elemental function clip_sp(x, xmin, xmax) result(res)
        real(sp), intent(in) :: x
        real(sp), intent(in) :: xmin
        real(sp), intent(in) :: xmax
        real(sp) :: res

        res = max(min(x, xmax), xmin)
    end function clip_sp

    elemental function clip_dp(x, xmin, xmax) result(res)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: xmin
        real(dp), intent(in) :: xmax
        real(dp) :: res

        res = max(min(x, xmax), xmin)
    end function clip_dp


    elemental function arg_sp(z) result(result) 
        complex(sp), intent(in) :: z
        real(sp) :: result

        result = merge(0.0_sp, atan2(z%im, z%re), z == (0.0_sp, 0.0_sp))

    end function arg_sp

    elemental function argd_sp(z) result(result) 
        complex(sp), intent(in) :: z
        real(sp) :: result

        result = merge(0.0_sp, atan2(z%im, z%re)*180.0_sp/PI_sp, &
                 z == (0.0_sp, 0.0_sp))

    end function argd_sp

    elemental function argpi_sp(z) result(result) 
        complex(sp), intent(in) :: z
        real(sp) :: result

        result = merge(0.0_sp, atan2(z%im, z%re)/PI_sp, &
                 z == (0.0_sp, 0.0_sp))
                 

    end function argpi_sp
    elemental function arg_dp(z) result(result) 
        complex(dp), intent(in) :: z
        real(dp) :: result

        result = merge(0.0_dp, atan2(z%im, z%re), z == (0.0_dp, 0.0_dp))

    end function arg_dp

    elemental function argd_dp(z) result(result) 
        complex(dp), intent(in) :: z
        real(dp) :: result

        result = merge(0.0_dp, atan2(z%im, z%re)*180.0_dp/PI_dp, &
                 z == (0.0_dp, 0.0_dp))

    end function argd_dp

    elemental function argpi_dp(z) result(result) 
        complex(dp), intent(in) :: z
        real(dp) :: result

        result = merge(0.0_dp, atan2(z%im, z%re)/PI_dp, &
                 z == (0.0_dp, 0.0_dp))
                 

    end function argpi_dp

    !> Returns the greatest common divisor of two integers of kind int8
    !> using the Euclidean algorithm.
    elemental function gcd_int8(a, b) result(res)
        integer(int8), intent(in) :: a
        integer(int8), intent(in) :: b
        integer(int8) :: res

        integer(int8) :: rem, tmp

        rem = min(abs(a), abs(b))
        res = max(abs(a), abs(b))
        do while (rem /= 0_int8)
          tmp = rem
          rem = mod(res, rem)
          res = tmp
        end do
    end function gcd_int8

    !> Returns the greatest common divisor of two integers of kind int16
    !> using the Euclidean algorithm.
    elemental function gcd_int16(a, b) result(res)
        integer(int16), intent(in) :: a
        integer(int16), intent(in) :: b
        integer(int16) :: res

        integer(int16) :: rem, tmp

        rem = min(abs(a), abs(b))
        res = max(abs(a), abs(b))
        do while (rem /= 0_int16)
          tmp = rem
          rem = mod(res, rem)
          res = tmp
        end do
    end function gcd_int16

    !> Returns the greatest common divisor of two integers of kind int32
    !> using the Euclidean algorithm.
    elemental function gcd_int32(a, b) result(res)
        integer(int32), intent(in) :: a
        integer(int32), intent(in) :: b
        integer(int32) :: res

        integer(int32) :: rem, tmp

        rem = min(abs(a), abs(b))
        res = max(abs(a), abs(b))
        do while (rem /= 0_int32)
          tmp = rem
          rem = mod(res, rem)
          res = tmp
        end do
    end function gcd_int32

    !> Returns the greatest common divisor of two integers of kind int64
    !> using the Euclidean algorithm.
    elemental function gcd_int64(a, b) result(res)
        integer(int64), intent(in) :: a
        integer(int64), intent(in) :: b
        integer(int64) :: res

        integer(int64) :: rem, tmp

        rem = min(abs(a), abs(b))
        res = max(abs(a), abs(b))
        do while (rem /= 0_int64)
          tmp = rem
          rem = mod(res, rem)
          res = tmp
        end do
    end function gcd_int64


    !> `arange` creates a vector of the `real(sp)` type 
    !>  with evenly spaced values within a given interval.
    pure module function arange_r_sp(start, end, step) result(result)

        real(sp), intent(in) :: start
        real(sp), intent(in), optional :: end, step
        real(sp), allocatable :: result(:)
        
        real(sp) :: start_, end_, step_
        integer :: i

        start_ = merge(start, 1.0_sp, present(end))
        end_   = optval(end, start)
        step_  = optval(step, 1.0_sp)
        step_  = sign(merge(step_, 1.0_sp, step_ /= 0.0_sp), end_ - start_)

        allocate(result(floor((end_ - start_)/step_) + 1))

        result = [(start_ + (i - 1)*step_, i=1, size(result), 1)]

    end function arange_r_sp
    !> `arange` creates a vector of the `real(dp)` type 
    !>  with evenly spaced values within a given interval.
    pure module function arange_r_dp(start, end, step) result(result)

        real(dp), intent(in) :: start
        real(dp), intent(in), optional :: end, step
        real(dp), allocatable :: result(:)
        
        real(dp) :: start_, end_, step_
        integer :: i

        start_ = merge(start, 1.0_dp, present(end))
        end_   = optval(end, start)
        step_  = optval(step, 1.0_dp)
        step_  = sign(merge(step_, 1.0_dp, step_ /= 0.0_dp), end_ - start_)

        allocate(result(floor((end_ - start_)/step_) + 1))

        result = [(start_ + (i - 1)*step_, i=1, size(result), 1)]

    end function arange_r_dp

    !> `arange` creates a vector of the `integer(int8)` type 
    !>  with evenly spaced values within a given interval.
    pure module function arange_i_int8(start, end, step) result(result)

        integer(int8), intent(in) :: start
        integer(int8), intent(in), optional :: end, step
        integer(int8), allocatable :: result(:)
        
        integer(int8) :: start_, end_, step_
        integer(int8) :: i

        start_ = merge(start, 1_int8, present(end))
        end_   = optval(end, start)
        step_  = optval(step, 1_int8)
        step_  = sign(merge(step_, 1_int8, step_ /= 0_int8), end_ - start_)

        allocate(result((end_ - start_)/step_ + 1_int8))

        result = [(i, i=start_, end_, step_)]

    end function arange_i_int8
    !> `arange` creates a vector of the `integer(int16)` type 
    !>  with evenly spaced values within a given interval.
    pure module function arange_i_int16(start, end, step) result(result)

        integer(int16), intent(in) :: start
        integer(int16), intent(in), optional :: end, step
        integer(int16), allocatable :: result(:)
        
        integer(int16) :: start_, end_, step_
        integer(int16) :: i

        start_ = merge(start, 1_int16, present(end))
        end_   = optval(end, start)
        step_  = optval(step, 1_int16)
        step_  = sign(merge(step_, 1_int16, step_ /= 0_int16), end_ - start_)

        allocate(result((end_ - start_)/step_ + 1_int16))

        result = [(i, i=start_, end_, step_)]

    end function arange_i_int16
    !> `arange` creates a vector of the `integer(int32)` type 
    !>  with evenly spaced values within a given interval.
    pure module function arange_i_int32(start, end, step) result(result)

        integer(int32), intent(in) :: start
        integer(int32), intent(in), optional :: end, step
        integer(int32), allocatable :: result(:)
        
        integer(int32) :: start_, end_, step_
        integer(int32) :: i

        start_ = merge(start, 1_int32, present(end))
        end_   = optval(end, start)
        step_  = optval(step, 1_int32)
        step_  = sign(merge(step_, 1_int32, step_ /= 0_int32), end_ - start_)

        allocate(result((end_ - start_)/step_ + 1_int32))

        result = [(i, i=start_, end_, step_)]

    end function arange_i_int32
    !> `arange` creates a vector of the `integer(int64)` type 
    !>  with evenly spaced values within a given interval.
    pure module function arange_i_int64(start, end, step) result(result)

        integer(int64), intent(in) :: start
        integer(int64), intent(in), optional :: end, step
        integer(int64), allocatable :: result(:)
        
        integer(int64) :: start_, end_, step_
        integer(int64) :: i

        start_ = merge(start, 1_int64, present(end))
        end_   = optval(end, start)
        step_  = optval(step, 1_int64)
        step_  = sign(merge(step_, 1_int64, step_ /= 0_int64), end_ - start_)

        allocate(result((end_ - start_)/step_ + 1_int64))

        result = [(i, i=start_, end_, step_)]

    end function arange_i_int64

    elemental module logical function is_close_rsp(a, b, rel_tol, abs_tol, equal_nan) result(close)
        use ieee_arithmetic, only: ieee_is_nan
        real(sp), intent(in) :: a, b
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan
        logical :: equal_nan_

        equal_nan_ = optval(equal_nan, .false.)
        
        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
        else
            close = abs(a - b) <= max( abs(optval(rel_tol, sqrt_eps_sp)*max(abs(a), abs(b))), &
                                       abs(optval(abs_tol, 0.0_sp)) )
        end if     

    end function is_close_rsp
    elemental module logical function is_close_rdp(a, b, rel_tol, abs_tol, equal_nan) result(close)
        use ieee_arithmetic, only: ieee_is_nan
        real(dp), intent(in) :: a, b
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan
        logical :: equal_nan_

        equal_nan_ = optval(equal_nan, .false.)
        
        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
        else
            close = abs(a - b) <= max( abs(optval(rel_tol, sqrt_eps_dp)*max(abs(a), abs(b))), &
                                       abs(optval(abs_tol, 0.0_dp)) )
        end if     

    end function is_close_rdp

    elemental module logical function is_close_csp(a, b, rel_tol, abs_tol, equal_nan) result(close)
        use ieee_arithmetic, only: ieee_is_nan
        complex(sp), intent(in) :: a, b
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = is_close_rsp(a%re, b%re, rel_tol, abs_tol, equal_nan) .and. &
                is_close_rsp(a%im, b%im, rel_tol, abs_tol, equal_nan)

    end function is_close_csp
    elemental module logical function is_close_cdp(a, b, rel_tol, abs_tol, equal_nan) result(close)
        use ieee_arithmetic, only: ieee_is_nan
        complex(dp), intent(in) :: a, b
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = is_close_rdp(a%re, b%re, rel_tol, abs_tol, equal_nan) .and. &
                is_close_rdp(a%im, b%im, rel_tol, abs_tol, equal_nan)

    end function is_close_cdp

    logical pure function all_close_1_rsp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(sp), intent(in) :: a(:), b(:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_1_rsp
    logical pure function all_close_2_rsp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(sp), intent(in) :: a(:,:), b(:,:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_2_rsp
    logical pure function all_close_3_rsp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(sp), intent(in) :: a(:,:,:), b(:,:,:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_3_rsp
    logical pure function all_close_4_rsp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(sp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_4_rsp
    logical pure function all_close_1_rdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(dp), intent(in) :: a(:), b(:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_1_rdp
    logical pure function all_close_2_rdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_2_rdp
    logical pure function all_close_3_rdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(dp), intent(in) :: a(:,:,:), b(:,:,:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_3_rdp
    logical pure function all_close_4_rdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        real(dp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_4_rdp
    logical pure function all_close_1_csp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(sp), intent(in) :: a(:), b(:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_1_csp
    logical pure function all_close_2_csp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(sp), intent(in) :: a(:,:), b(:,:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_2_csp
    logical pure function all_close_3_csp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(sp), intent(in) :: a(:,:,:), b(:,:,:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_3_csp
    logical pure function all_close_4_csp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(sp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(sp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_4_csp
    logical pure function all_close_1_cdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(dp), intent(in) :: a(:), b(:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_1_cdp
    logical pure function all_close_2_cdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_2_cdp
    logical pure function all_close_3_cdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(dp), intent(in) :: a(:,:,:), b(:,:,:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_3_cdp
    logical pure function all_close_4_cdp(a, b, rel_tol, abs_tol, equal_nan) result(close)

        complex(dp), intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(dp), intent(in), optional :: rel_tol, abs_tol
        logical, intent(in), optional :: equal_nan

        close = all(is_close(a, b, rel_tol, abs_tol, equal_nan))

    end function all_close_4_cdp

    pure module function linspace_default_1_rsp_rsp(start, end) result(res)
      real(sp), intent(in) :: start
      real(sp), intent(in) :: end

      real(sp) :: res(DEFAULT_LINSPACE_LENGTH)

      res = linspace(start, end, DEFAULT_LINSPACE_LENGTH)

    end function linspace_default_1_rsp_rsp
    pure module function linspace_default_1_rdp_rdp(start, end) result(res)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end

      real(dp) :: res(DEFAULT_LINSPACE_LENGTH)

      res = linspace(start, end, DEFAULT_LINSPACE_LENGTH)

    end function linspace_default_1_rdp_rdp

    pure module function linspace_n_1_rsp_rsp(start, end_, n) result(res)
      real(sp), intent(in) :: start
      real(sp), intent(in) :: end_
      integer, intent(in) :: n

      real(sp) :: res(max(n, 0))

      integer :: i    ! Looping index
      real(sp) :: interval ! Difference between adjacent elements


      if(n <= 0) return ! If passed length is less than or equal to 0, return an empty (allocated with length 0) array
      if(n == 1) then
        res(1) = end_
        return
      end if

      interval = (end_ - start) / real((n - 1), sp)

      res(1) = start
      res(n) = end_

      do i = 2, n - 1

        res(i) = real((i-1), sp) * interval + start

      end do

    end function linspace_n_1_rsp_rsp
    pure module function linspace_n_1_rdp_rdp(start, end_, n) result(res)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end_
      integer, intent(in) :: n

      real(dp) :: res(max(n, 0))

      integer :: i    ! Looping index
      real(dp) :: interval ! Difference between adjacent elements


      if(n <= 0) return ! If passed length is less than or equal to 0, return an empty (allocated with length 0) array
      if(n == 1) then
        res(1) = end_
        return
      end if

      interval = (end_ - start) / real((n - 1), dp)

      res(1) = start
      res(n) = end_

      do i = 2, n - 1

        res(i) = real((i-1), dp) * interval + start

      end do

    end function linspace_n_1_rdp_rdp

      pure module function linspace_n_1_csp_csp(start, end, n) result(res)
        complex(sp), intent(in) :: start
        complex(sp), intent(in) :: end
        integer, intent(in) :: n

        complex(sp) :: res(max(n, 0))

        real(sp) :: x(max(n, 0)) ! array of the real part of complex number
        real(sp) :: y(max(n, 0)) ! array of the imaginary part of the complex number

        x = linspace(start%re, end%re, n)
        y = linspace(start%im, end%im, n)

        res = cmplx(x, y, kind=sp)
      end function linspace_n_1_csp_csp
      pure module function linspace_n_1_cdp_cdp(start, end, n) result(res)
        complex(dp), intent(in) :: start
        complex(dp), intent(in) :: end
        integer, intent(in) :: n

        complex(dp) :: res(max(n, 0))

        real(dp) :: x(max(n, 0)) ! array of the real part of complex number
        real(dp) :: y(max(n, 0)) ! array of the imaginary part of the complex number

        x = linspace(start%re, end%re, n)
        y = linspace(start%im, end%im, n)

        res = cmplx(x, y, kind=dp)
      end function linspace_n_1_cdp_cdp

    ! Add support for integer linspace
    !!
    !! When dealing with integers as the `start` and `end` parameters, the return type is always a `real(dp)`.
      pure module function linspace_default_1_iint8_iint8(start, end) result(res)
        integer(int8), intent(in) :: start
        integer(int8), intent(in) :: end

        real(dp) :: res(DEFAULT_LINSPACE_LENGTH)

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)
      end function linspace_default_1_iint8_iint8
      pure module function linspace_default_1_iint16_iint16(start, end) result(res)
        integer(int16), intent(in) :: start
        integer(int16), intent(in) :: end

        real(dp) :: res(DEFAULT_LINSPACE_LENGTH)

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)
      end function linspace_default_1_iint16_iint16
      pure module function linspace_default_1_iint32_iint32(start, end) result(res)
        integer(int32), intent(in) :: start
        integer(int32), intent(in) :: end

        real(dp) :: res(DEFAULT_LINSPACE_LENGTH)

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)
      end function linspace_default_1_iint32_iint32
      pure module function linspace_default_1_iint64_iint64(start, end) result(res)
        integer(int64), intent(in) :: start
        integer(int64), intent(in) :: end

        real(dp) :: res(DEFAULT_LINSPACE_LENGTH)

        res = linspace(real(start, kind=dp), real(end, kind=dp), DEFAULT_LINSPACE_LENGTH)
      end function linspace_default_1_iint64_iint64

      pure module function linspace_n_1_iint8_iint8(start, end, n) result(res)
        integer(int8), intent(in) :: start
        integer(int8), intent(in) :: end
        integer, intent(in) :: n

        real(dp) :: res(max(n, 0))

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end function linspace_n_1_iint8_iint8
      pure module function linspace_n_1_iint16_iint16(start, end, n) result(res)
        integer(int16), intent(in) :: start
        integer(int16), intent(in) :: end
        integer, intent(in) :: n

        real(dp) :: res(max(n, 0))

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end function linspace_n_1_iint16_iint16
      pure module function linspace_n_1_iint32_iint32(start, end, n) result(res)
        integer(int32), intent(in) :: start
        integer(int32), intent(in) :: end
        integer, intent(in) :: n

        real(dp) :: res(max(n, 0))

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end function linspace_n_1_iint32_iint32
      pure module function linspace_n_1_iint64_iint64(start, end, n) result(res)
        integer(int64), intent(in) :: start
        integer(int64), intent(in) :: end
        integer, intent(in) :: n

        real(dp) :: res(max(n, 0))

        res = linspace(real(start, kind=dp), real(end, kind=dp), n)

      end function linspace_n_1_iint64_iint64

      pure module function linspace_default_1_csp_csp(start, end) result(res)
        complex(sp), intent(in) :: start
        complex(sp), intent(in) :: end

        complex(sp) :: res(DEFAULT_LINSPACE_LENGTH)

        res = linspace(start, end, DEFAULT_LINSPACE_LENGTH)
      end function linspace_default_1_csp_csp
      pure module function linspace_default_1_cdp_cdp(start, end) result(res)
        complex(dp), intent(in) :: start
        complex(dp), intent(in) :: end

        complex(dp) :: res(DEFAULT_LINSPACE_LENGTH)

        res = linspace(start, end, DEFAULT_LINSPACE_LENGTH)
      end function linspace_default_1_cdp_cdp

    pure module function logspace_1_rsp_default(start, end) result(res)

      real(sp), intent(in) :: start
      real(sp), intent(in) :: end

      real(sp) :: res(DEFAULT_LOGSPACE_LENGTH)

      res = logspace(start, end, DEFAULT_LOGSPACE_LENGTH, real(DEFAULT_LOGSPACE_BASE, sp))

    end function logspace_1_rsp_default
    pure module function logspace_1_rdp_default(start, end) result(res)

      real(dp), intent(in) :: start
      real(dp), intent(in) :: end

      real(dp) :: res(DEFAULT_LOGSPACE_LENGTH)

      res = logspace(start, end, DEFAULT_LOGSPACE_LENGTH, real(DEFAULT_LOGSPACE_BASE, dp))

    end function logspace_1_rdp_default
    pure module function logspace_1_csp_default(start, end) result(res)

      complex(sp), intent(in) :: start
      complex(sp), intent(in) :: end

      complex(sp) :: res(DEFAULT_LOGSPACE_LENGTH)

      res = logspace(start, end, DEFAULT_LOGSPACE_LENGTH, real(DEFAULT_LOGSPACE_BASE, sp))

    end function logspace_1_csp_default
    pure module function logspace_1_cdp_default(start, end) result(res)

      complex(dp), intent(in) :: start
      complex(dp), intent(in) :: end

      complex(dp) :: res(DEFAULT_LOGSPACE_LENGTH)

      res = logspace(start, end, DEFAULT_LOGSPACE_LENGTH, real(DEFAULT_LOGSPACE_BASE, dp))

    end function logspace_1_cdp_default

    pure module function logspace_1_iint32_default(start, end) result(res)

      integer, intent(in) :: start
      integer, intent(in) :: end

      real(dp) :: res(DEFAULT_LOGSPACE_LENGTH)

      res = logspace(start, end, DEFAULT_LOGSPACE_LENGTH, DEFAULT_LOGSPACE_BASE)

  end function logspace_1_iint32_default

    pure module function logspace_1_rsp_n(start, end, n) result(res)
      real(sp), intent(in) :: start
      real(sp), intent(in) :: end
      integer, intent(in) :: n

      real(sp) :: res(max(n, 0))

      res = logspace(start, end, n, real(DEFAULT_LOGSPACE_BASE, sp))
    end function logspace_1_rsp_n
    pure module function logspace_1_rdp_n(start, end, n) result(res)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end
      integer, intent(in) :: n

      real(dp) :: res(max(n, 0))

      res = logspace(start, end, n, real(DEFAULT_LOGSPACE_BASE, dp))
    end function logspace_1_rdp_n
    pure module function logspace_1_csp_n(start, end, n) result(res)
      complex(sp), intent(in) :: start
      complex(sp), intent(in) :: end
      integer, intent(in) :: n

      complex(sp) :: res(max(n, 0))

      res = logspace(start, end, n, real(DEFAULT_LOGSPACE_BASE, sp))
    end function logspace_1_csp_n
    pure module function logspace_1_cdp_n(start, end, n) result(res)
      complex(dp), intent(in) :: start
      complex(dp), intent(in) :: end
      integer, intent(in) :: n

      complex(dp) :: res(max(n, 0))

      res = logspace(start, end, n, real(DEFAULT_LOGSPACE_BASE, dp))
    end function logspace_1_cdp_n

    pure module function logspace_1_iint32_n(start, end, n) result(res)
      integer, intent(in) :: start
      integer, intent(in) :: end
      integer, intent(in) :: n

      real(dp) :: res(n)

      res = logspace(start, end, n, DEFAULT_LOGSPACE_BASE)
    end function logspace_1_iint32_n

    ! Generate logarithmically spaced sequence from sp base to the powers
    ! of sp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    pure module function logspace_1_rsp_n_rbase(start, end, n, base) result(res)
      real(sp), intent(in) :: start
      real(sp), intent(in) :: end
      integer, intent(in) :: n
      real(sp), intent(in) :: base
      ! real(sp) endpoints + real(sp) base = real(sp) result
      real(sp) :: res(max(n, 0))

      real(sp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_rsp_n_rbase

    pure module function logspace_1_rsp_n_cbase(start, end, n, base) result(res)
      real(sp), intent(in) :: start
      real(sp), intent(in) :: end
      integer, intent(in) :: n
      complex(sp), intent(in) :: base
      ! real(sp) endpoints + complex(sp) base = complex(sp) result
      real(sp) :: res(max(n, 0))

      real(sp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_rsp_n_cbase

    pure module function logspace_1_rsp_n_ibase(start, end, n, base) result(res)
      real(sp), intent(in) :: start
      real(sp), intent(in) :: end
      integer, intent(in) :: n
      integer, intent(in) :: base
      ! real(sp) endpoints + integer base = real(sp) result
      real(sp) :: res(max(n, 0))

      real(sp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_rsp_n_ibase
    ! Generate logarithmically spaced sequence from dp base to the powers
    ! of dp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    pure module function logspace_1_rdp_n_rbase(start, end, n, base) result(res)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end
      integer, intent(in) :: n
      real(dp), intent(in) :: base
      ! real(dp) endpoints + real(dp) base = real(dp) result
      real(dp) :: res(max(n, 0))

      real(dp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_rdp_n_rbase

    pure module function logspace_1_rdp_n_cbase(start, end, n, base) result(res)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end
      integer, intent(in) :: n
      complex(dp), intent(in) :: base
      ! real(dp) endpoints + complex(dp) base = complex(dp) result
      real(dp) :: res(max(n, 0))

      real(dp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_rdp_n_cbase

    pure module function logspace_1_rdp_n_ibase(start, end, n, base) result(res)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end
      integer, intent(in) :: n
      integer, intent(in) :: base
      ! real(dp) endpoints + integer base = real(dp) result
      real(dp) :: res(max(n, 0))

      real(dp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_rdp_n_ibase

    ! Generate logarithmically spaced sequence from sp base to the powers
    ! of sp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    pure module function logspace_1_csp_n_rbase(start, end, n, base) result(res)
      complex(sp), intent(in) :: start
      complex(sp), intent(in) :: end
      integer, intent(in) :: n
      real(sp), intent(in) :: base
      ! complex(sp) endpoints + real(sp) base = complex(sp) result
      complex(sp) :: res(max(n, 0))

      complex(sp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_csp_n_rbase

    pure module function logspace_1_csp_n_cbase(start, end, n, base) result(res)
      complex(sp), intent(in) :: start
      complex(sp), intent(in) :: end
      integer, intent(in) :: n
      complex(sp), intent(in) :: base
      ! complex(sp) endpoints + complex(sp) base = complex(sp) result
      complex(sp) :: res(max(n, 0))

      complex(sp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_csp_n_cbase

    pure module function logspace_1_csp_n_ibase(start, end, n, base) result(res)
      complex(sp), intent(in) :: start
      complex(sp), intent(in) :: end
      integer, intent(in) :: n
      integer, intent(in) :: base
      ! complex(sp) endpoints + integer base = complex(sp) result
      complex(sp) :: res(max(n, 0))

      complex(sp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_csp_n_ibase
    ! Generate logarithmically spaced sequence from dp base to the powers
    ! of dp start and end. [base^start, ... , base^end]
    ! Different combinations of parameter types will lead to different result types.
    ! Those combinations are indicated in the body of each function.
    pure module function logspace_1_cdp_n_rbase(start, end, n, base) result(res)
      complex(dp), intent(in) :: start
      complex(dp), intent(in) :: end
      integer, intent(in) :: n
      real(dp), intent(in) :: base
      ! complex(dp) endpoints + real(dp) base = complex(dp) result
      complex(dp) :: res(max(n, 0))

      complex(dp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_cdp_n_rbase

    pure module function logspace_1_cdp_n_cbase(start, end, n, base) result(res)
      complex(dp), intent(in) :: start
      complex(dp), intent(in) :: end
      integer, intent(in) :: n
      complex(dp), intent(in) :: base
      ! complex(dp) endpoints + complex(dp) base = complex(dp) result
      complex(dp) :: res(max(n, 0))

      complex(dp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_cdp_n_cbase

    pure module function logspace_1_cdp_n_ibase(start, end, n, base) result(res)
      complex(dp), intent(in) :: start
      complex(dp), intent(in) :: end
      integer, intent(in) :: n
      integer, intent(in) :: base
      ! complex(dp) endpoints + integer base = complex(dp) result
      complex(dp) :: res(max(n, 0))

      complex(dp) :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_cdp_n_ibase

    pure module function logspace_1_iint32_n_rspbase(start, end, n, base) result(res)
      integer, intent(in) :: start
      integer, intent(in) :: end
      integer, intent(in) :: n
      real(sp), intent(in) :: base
      ! integer endpoints + real(sp) base = real(sp) result
      real(sp) :: res(max(n, 0))

      integer :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_iint32_n_rspbase

    pure module function logspace_1_iint32_n_cspbase(start, end, n, base) result(res)
      integer, intent(in) :: start
      integer, intent(in) :: end
      integer, intent(in) :: n
      complex(sp), intent(in) :: base
      ! integer endpoints + complex(sp) base = complex(sp) result
      complex(sp) :: res(max(n, 0))

      integer :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_iint32_n_cspbase
    pure module function logspace_1_iint32_n_rdpbase(start, end, n, base) result(res)
      integer, intent(in) :: start
      integer, intent(in) :: end
      integer, intent(in) :: n
      real(dp), intent(in) :: base
      ! integer endpoints + real(dp) base = real(dp) result
      real(dp) :: res(max(n, 0))

      integer :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_iint32_n_rdpbase

    pure module function logspace_1_iint32_n_cdpbase(start, end, n, base) result(res)
      integer, intent(in) :: start
      integer, intent(in) :: end
      integer, intent(in) :: n
      complex(dp), intent(in) :: base
      ! integer endpoints + complex(dp) base = complex(dp) result
      complex(dp) :: res(max(n, 0))

      integer :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_iint32_n_cdpbase

    pure module function logspace_1_iint32_n_ibase(start, end, n, base) result(res)
      integer, intent(in) :: start
      integer, intent(in) :: end
      integer, intent(in) :: n
      integer, intent(in) :: base
      ! integer endpoints + integer base = integer result
      integer :: res(max(n, 0))

      integer :: exponents(max(n, 0))
      exponents = linspace(start, end, n)
      res = base ** exponents
    end function logspace_1_iint32_n_ibase

    !> `diff` computes differences of adjacent elements of an array.
    
    pure module function diff_1_sp(x, n, prepend, append) result(y)
        real(sp), intent(in) :: x(:)
        integer, intent(in), optional :: n
        real(sp), intent(in), optional :: prepend(:), append(:)
        real(sp), allocatable :: y(:)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(prepend)) size_prepend = size(prepend) 
        if (present(append)) size_append = size(append)
        size_x = size(x)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0))
            return
        end if

        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size_x-1))
            y = x(2:) - x(1:size_x-1)
            return
        end if

        block
        real(sp) :: work(size_work)
        if (size_prepend > 0) work(:size_prepend) = prepend
        work(size_prepend+1:size_prepend+size_x) = x
        if (size_append > 0) work(size_prepend+size_x+1:) = append
        
        do i = 1, n_
            work(1:size_work-i) = work(2:size_work-i+1) - work(1:size_work-i)
        end do
        allocate(y(size_work-n_))
        y = work(1:size_work-n_)
        end block

    end function diff_1_sp

    pure module function diff_2_sp(x, n, dim, prepend, append) result(y)
        real(sp), intent(in) :: x(:, :)
        integer, intent(in), optional :: n, dim
        real(sp), intent(in), optional :: prepend(:, :), append(:, :)
        real(sp), allocatable :: y(:, :)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, dim_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(dim)) then
            if (dim == 1 .or. dim == 2) then
                dim_ = dim
            else
                dim_ = 1
            end if
        else
            dim_ = 1
        end if
        
        if (present(prepend)) size_prepend = size(prepend, dim_)
        if (present(append)) size_append = size(append, dim_)
        size_x = size(x, dim_)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0, 0))
            return
        end if
        
        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size(x, 1), size(x, 2)))
            if (dim_ == 1) then
                y = x(2:, :) - x(1:size_x-1, :)
            elseif (dim_ == 2) then
                y = x(:, 2:) - x(:, 1:size_x-1)
            end if
            return
        end if
        
        if (dim_ == 1) then
            block
            real(sp) :: work(size_work, size(x, 2))
            if (size_prepend > 0) work(1:size_prepend, :) = prepend
            work(size_prepend+1:size_x+size_prepend, :) = x
            if (size_append > 0) work(size_x+size_prepend+1:, :) = append
            do i = 1, n_
                work(1:size_work-i, :) = work(2:size_work-i+1, :) - work(1:size_work-i, :)
            end do
            allocate(y(size_work-n_, size(x, 2)))
            y = work(1:size_work-n_, :)
            end block
            
        elseif (dim_ == 2) then
            block
            real(sp) :: work(size(x, 1), size_work)
            if (size_prepend > 0) work(:, 1:size_prepend) = prepend
            work(:, size_prepend+1:size_x+size_prepend) = x
            if (size_append > 0) work(:, size_x+size_prepend+1:) = append
            do i = 1, n_
                work(:, 1:size_work-i) = work(:, 2:size_work-i+1) - work(:, 1:size_work-i)
            end do
            allocate(y(size(x, 1), size_work-n_))
            y = work(:, 1:size_work-n_)
            end block
            
        end if

    end function diff_2_sp
    pure module function diff_1_dp(x, n, prepend, append) result(y)
        real(dp), intent(in) :: x(:)
        integer, intent(in), optional :: n
        real(dp), intent(in), optional :: prepend(:), append(:)
        real(dp), allocatable :: y(:)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(prepend)) size_prepend = size(prepend) 
        if (present(append)) size_append = size(append)
        size_x = size(x)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0))
            return
        end if

        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size_x-1))
            y = x(2:) - x(1:size_x-1)
            return
        end if

        block
        real(dp) :: work(size_work)
        if (size_prepend > 0) work(:size_prepend) = prepend
        work(size_prepend+1:size_prepend+size_x) = x
        if (size_append > 0) work(size_prepend+size_x+1:) = append
        
        do i = 1, n_
            work(1:size_work-i) = work(2:size_work-i+1) - work(1:size_work-i)
        end do
        allocate(y(size_work-n_))
        y = work(1:size_work-n_)
        end block

    end function diff_1_dp

    pure module function diff_2_dp(x, n, dim, prepend, append) result(y)
        real(dp), intent(in) :: x(:, :)
        integer, intent(in), optional :: n, dim
        real(dp), intent(in), optional :: prepend(:, :), append(:, :)
        real(dp), allocatable :: y(:, :)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, dim_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(dim)) then
            if (dim == 1 .or. dim == 2) then
                dim_ = dim
            else
                dim_ = 1
            end if
        else
            dim_ = 1
        end if
        
        if (present(prepend)) size_prepend = size(prepend, dim_)
        if (present(append)) size_append = size(append, dim_)
        size_x = size(x, dim_)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0, 0))
            return
        end if
        
        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size(x, 1), size(x, 2)))
            if (dim_ == 1) then
                y = x(2:, :) - x(1:size_x-1, :)
            elseif (dim_ == 2) then
                y = x(:, 2:) - x(:, 1:size_x-1)
            end if
            return
        end if
        
        if (dim_ == 1) then
            block
            real(dp) :: work(size_work, size(x, 2))
            if (size_prepend > 0) work(1:size_prepend, :) = prepend
            work(size_prepend+1:size_x+size_prepend, :) = x
            if (size_append > 0) work(size_x+size_prepend+1:, :) = append
            do i = 1, n_
                work(1:size_work-i, :) = work(2:size_work-i+1, :) - work(1:size_work-i, :)
            end do
            allocate(y(size_work-n_, size(x, 2)))
            y = work(1:size_work-n_, :)
            end block
            
        elseif (dim_ == 2) then
            block
            real(dp) :: work(size(x, 1), size_work)
            if (size_prepend > 0) work(:, 1:size_prepend) = prepend
            work(:, size_prepend+1:size_x+size_prepend) = x
            if (size_append > 0) work(:, size_x+size_prepend+1:) = append
            do i = 1, n_
                work(:, 1:size_work-i) = work(:, 2:size_work-i+1) - work(:, 1:size_work-i)
            end do
            allocate(y(size(x, 1), size_work-n_))
            y = work(:, 1:size_work-n_)
            end block
            
        end if

    end function diff_2_dp
    pure module function diff_1_int8(x, n, prepend, append) result(y)
        integer(int8), intent(in) :: x(:)
        integer, intent(in), optional :: n
        integer(int8), intent(in), optional :: prepend(:), append(:)
        integer(int8), allocatable :: y(:)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(prepend)) size_prepend = size(prepend) 
        if (present(append)) size_append = size(append)
        size_x = size(x)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0))
            return
        end if

        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size_x-1))
            y = x(2:) - x(1:size_x-1)
            return
        end if

        block
        integer(int8) :: work(size_work)
        if (size_prepend > 0) work(:size_prepend) = prepend
        work(size_prepend+1:size_prepend+size_x) = x
        if (size_append > 0) work(size_prepend+size_x+1:) = append
        
        do i = 1, n_
            work(1:size_work-i) = work(2:size_work-i+1) - work(1:size_work-i)
        end do
        allocate(y(size_work-n_))
        y = work(1:size_work-n_)
        end block

    end function diff_1_int8

    pure module function diff_2_int8(x, n, dim, prepend, append) result(y)
        integer(int8), intent(in) :: x(:, :)
        integer, intent(in), optional :: n, dim
        integer(int8), intent(in), optional :: prepend(:, :), append(:, :)
        integer(int8), allocatable :: y(:, :)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, dim_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(dim)) then
            if (dim == 1 .or. dim == 2) then
                dim_ = dim
            else
                dim_ = 1
            end if
        else
            dim_ = 1
        end if
        
        if (present(prepend)) size_prepend = size(prepend, dim_)
        if (present(append)) size_append = size(append, dim_)
        size_x = size(x, dim_)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0, 0))
            return
        end if
        
        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size(x, 1), size(x, 2)))
            if (dim_ == 1) then
                y = x(2:, :) - x(1:size_x-1, :)
            elseif (dim_ == 2) then
                y = x(:, 2:) - x(:, 1:size_x-1)
            end if
            return
        end if
        
        if (dim_ == 1) then
            block
            integer(int8) :: work(size_work, size(x, 2))
            if (size_prepend > 0) work(1:size_prepend, :) = prepend
            work(size_prepend+1:size_x+size_prepend, :) = x
            if (size_append > 0) work(size_x+size_prepend+1:, :) = append
            do i = 1, n_
                work(1:size_work-i, :) = work(2:size_work-i+1, :) - work(1:size_work-i, :)
            end do
            allocate(y(size_work-n_, size(x, 2)))
            y = work(1:size_work-n_, :)
            end block
            
        elseif (dim_ == 2) then
            block
            integer(int8) :: work(size(x, 1), size_work)
            if (size_prepend > 0) work(:, 1:size_prepend) = prepend
            work(:, size_prepend+1:size_x+size_prepend) = x
            if (size_append > 0) work(:, size_x+size_prepend+1:) = append
            do i = 1, n_
                work(:, 1:size_work-i) = work(:, 2:size_work-i+1) - work(:, 1:size_work-i)
            end do
            allocate(y(size(x, 1), size_work-n_))
            y = work(:, 1:size_work-n_)
            end block
            
        end if

    end function diff_2_int8
    pure module function diff_1_int16(x, n, prepend, append) result(y)
        integer(int16), intent(in) :: x(:)
        integer, intent(in), optional :: n
        integer(int16), intent(in), optional :: prepend(:), append(:)
        integer(int16), allocatable :: y(:)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(prepend)) size_prepend = size(prepend) 
        if (present(append)) size_append = size(append)
        size_x = size(x)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0))
            return
        end if

        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size_x-1))
            y = x(2:) - x(1:size_x-1)
            return
        end if

        block
        integer(int16) :: work(size_work)
        if (size_prepend > 0) work(:size_prepend) = prepend
        work(size_prepend+1:size_prepend+size_x) = x
        if (size_append > 0) work(size_prepend+size_x+1:) = append
        
        do i = 1, n_
            work(1:size_work-i) = work(2:size_work-i+1) - work(1:size_work-i)
        end do
        allocate(y(size_work-n_))
        y = work(1:size_work-n_)
        end block

    end function diff_1_int16

    pure module function diff_2_int16(x, n, dim, prepend, append) result(y)
        integer(int16), intent(in) :: x(:, :)
        integer, intent(in), optional :: n, dim
        integer(int16), intent(in), optional :: prepend(:, :), append(:, :)
        integer(int16), allocatable :: y(:, :)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, dim_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(dim)) then
            if (dim == 1 .or. dim == 2) then
                dim_ = dim
            else
                dim_ = 1
            end if
        else
            dim_ = 1
        end if
        
        if (present(prepend)) size_prepend = size(prepend, dim_)
        if (present(append)) size_append = size(append, dim_)
        size_x = size(x, dim_)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0, 0))
            return
        end if
        
        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size(x, 1), size(x, 2)))
            if (dim_ == 1) then
                y = x(2:, :) - x(1:size_x-1, :)
            elseif (dim_ == 2) then
                y = x(:, 2:) - x(:, 1:size_x-1)
            end if
            return
        end if
        
        if (dim_ == 1) then
            block
            integer(int16) :: work(size_work, size(x, 2))
            if (size_prepend > 0) work(1:size_prepend, :) = prepend
            work(size_prepend+1:size_x+size_prepend, :) = x
            if (size_append > 0) work(size_x+size_prepend+1:, :) = append
            do i = 1, n_
                work(1:size_work-i, :) = work(2:size_work-i+1, :) - work(1:size_work-i, :)
            end do
            allocate(y(size_work-n_, size(x, 2)))
            y = work(1:size_work-n_, :)
            end block
            
        elseif (dim_ == 2) then
            block
            integer(int16) :: work(size(x, 1), size_work)
            if (size_prepend > 0) work(:, 1:size_prepend) = prepend
            work(:, size_prepend+1:size_x+size_prepend) = x
            if (size_append > 0) work(:, size_x+size_prepend+1:) = append
            do i = 1, n_
                work(:, 1:size_work-i) = work(:, 2:size_work-i+1) - work(:, 1:size_work-i)
            end do
            allocate(y(size(x, 1), size_work-n_))
            y = work(:, 1:size_work-n_)
            end block
            
        end if

    end function diff_2_int16
    pure module function diff_1_int32(x, n, prepend, append) result(y)
        integer(int32), intent(in) :: x(:)
        integer, intent(in), optional :: n
        integer(int32), intent(in), optional :: prepend(:), append(:)
        integer(int32), allocatable :: y(:)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(prepend)) size_prepend = size(prepend) 
        if (present(append)) size_append = size(append)
        size_x = size(x)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0))
            return
        end if

        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size_x-1))
            y = x(2:) - x(1:size_x-1)
            return
        end if

        block
        integer(int32) :: work(size_work)
        if (size_prepend > 0) work(:size_prepend) = prepend
        work(size_prepend+1:size_prepend+size_x) = x
        if (size_append > 0) work(size_prepend+size_x+1:) = append
        
        do i = 1, n_
            work(1:size_work-i) = work(2:size_work-i+1) - work(1:size_work-i)
        end do
        allocate(y(size_work-n_))
        y = work(1:size_work-n_)
        end block

    end function diff_1_int32

    pure module function diff_2_int32(x, n, dim, prepend, append) result(y)
        integer(int32), intent(in) :: x(:, :)
        integer, intent(in), optional :: n, dim
        integer(int32), intent(in), optional :: prepend(:, :), append(:, :)
        integer(int32), allocatable :: y(:, :)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, dim_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(dim)) then
            if (dim == 1 .or. dim == 2) then
                dim_ = dim
            else
                dim_ = 1
            end if
        else
            dim_ = 1
        end if
        
        if (present(prepend)) size_prepend = size(prepend, dim_)
        if (present(append)) size_append = size(append, dim_)
        size_x = size(x, dim_)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0, 0))
            return
        end if
        
        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size(x, 1), size(x, 2)))
            if (dim_ == 1) then
                y = x(2:, :) - x(1:size_x-1, :)
            elseif (dim_ == 2) then
                y = x(:, 2:) - x(:, 1:size_x-1)
            end if
            return
        end if
        
        if (dim_ == 1) then
            block
            integer(int32) :: work(size_work, size(x, 2))
            if (size_prepend > 0) work(1:size_prepend, :) = prepend
            work(size_prepend+1:size_x+size_prepend, :) = x
            if (size_append > 0) work(size_x+size_prepend+1:, :) = append
            do i = 1, n_
                work(1:size_work-i, :) = work(2:size_work-i+1, :) - work(1:size_work-i, :)
            end do
            allocate(y(size_work-n_, size(x, 2)))
            y = work(1:size_work-n_, :)
            end block
            
        elseif (dim_ == 2) then
            block
            integer(int32) :: work(size(x, 1), size_work)
            if (size_prepend > 0) work(:, 1:size_prepend) = prepend
            work(:, size_prepend+1:size_x+size_prepend) = x
            if (size_append > 0) work(:, size_x+size_prepend+1:) = append
            do i = 1, n_
                work(:, 1:size_work-i) = work(:, 2:size_work-i+1) - work(:, 1:size_work-i)
            end do
            allocate(y(size(x, 1), size_work-n_))
            y = work(:, 1:size_work-n_)
            end block
            
        end if

    end function diff_2_int32
    pure module function diff_1_int64(x, n, prepend, append) result(y)
        integer(int64), intent(in) :: x(:)
        integer, intent(in), optional :: n
        integer(int64), intent(in), optional :: prepend(:), append(:)
        integer(int64), allocatable :: y(:)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(prepend)) size_prepend = size(prepend) 
        if (present(append)) size_append = size(append)
        size_x = size(x)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0))
            return
        end if

        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size_x-1))
            y = x(2:) - x(1:size_x-1)
            return
        end if

        block
        integer(int64) :: work(size_work)
        if (size_prepend > 0) work(:size_prepend) = prepend
        work(size_prepend+1:size_prepend+size_x) = x
        if (size_append > 0) work(size_prepend+size_x+1:) = append
        
        do i = 1, n_
            work(1:size_work-i) = work(2:size_work-i+1) - work(1:size_work-i)
        end do
        allocate(y(size_work-n_))
        y = work(1:size_work-n_)
        end block

    end function diff_1_int64

    pure module function diff_2_int64(x, n, dim, prepend, append) result(y)
        integer(int64), intent(in) :: x(:, :)
        integer, intent(in), optional :: n, dim
        integer(int64), intent(in), optional :: prepend(:, :), append(:, :)
        integer(int64), allocatable :: y(:, :)
        integer :: size_prepend, size_append, size_x, size_work
        integer :: n_, dim_, i

        n_ = optval(n, 1)
        if (n_ <= 0) then
            y = x
            return
        end if
        
        size_prepend = 0
        size_append = 0
        if (present(dim)) then
            if (dim == 1 .or. dim == 2) then
                dim_ = dim
            else
                dim_ = 1
            end if
        else
            dim_ = 1
        end if
        
        if (present(prepend)) size_prepend = size(prepend, dim_)
        if (present(append)) size_append = size(append, dim_)
        size_x = size(x, dim_)
        size_work = size_x + size_prepend + size_append
        
        if (size_work <= n_) then
            allocate(y(0, 0))
            return
        end if
        
        !> Use a quick exit for the common case, to avoid memory allocation.
        if (size_prepend == 0 .and. size_append == 0 .and. n_ == 1) then
            allocate(y(size(x, 1), size(x, 2)))
            if (dim_ == 1) then
                y = x(2:, :) - x(1:size_x-1, :)
            elseif (dim_ == 2) then
                y = x(:, 2:) - x(:, 1:size_x-1)
            end if
            return
        end if
        
        if (dim_ == 1) then
            block
            integer(int64) :: work(size_work, size(x, 2))
            if (size_prepend > 0) work(1:size_prepend, :) = prepend
            work(size_prepend+1:size_x+size_prepend, :) = x
            if (size_append > 0) work(size_x+size_prepend+1:, :) = append
            do i = 1, n_
                work(1:size_work-i, :) = work(2:size_work-i+1, :) - work(1:size_work-i, :)
            end do
            allocate(y(size_work-n_, size(x, 2)))
            y = work(1:size_work-n_, :)
            end block
            
        elseif (dim_ == 2) then
            block
            integer(int64) :: work(size(x, 1), size_work)
            if (size_prepend > 0) work(:, 1:size_prepend) = prepend
            work(:, size_prepend+1:size_x+size_prepend) = x
            if (size_append > 0) work(:, size_x+size_prepend+1:) = append
            do i = 1, n_
                work(:, 1:size_work-i) = work(:, 2:size_work-i+1) - work(:, 1:size_work-i)
            end do
            allocate(y(size(x, 1), size_work-n_))
            y = work(:, 1:size_work-n_)
            end block
            
        end if

    end function diff_2_int64
end module stdlib_math
