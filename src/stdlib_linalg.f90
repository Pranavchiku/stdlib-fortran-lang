module stdlib_linalg
  !!Provides a support for various linear algebra procedures
  !! ([Specification](../page/specs/stdlib_linalg.html))
  use stdlib_kinds, only: sp, dp, xdp, qp, &
    int8, int16, int32, int64
  use stdlib_error, only: error_stop
  use stdlib_optval, only: optval
  implicit none
  private

  public :: diag
  public :: eye
  public :: trace
  public :: outer_product
  public :: kronecker_product
  public :: cross_product
  public :: is_square
  public :: is_diagonal
  public :: is_symmetric
  public :: is_skew_symmetric
  public :: is_hermitian
  public :: is_triangular
  public :: is_hessenberg

  interface diag
    !! version: experimental
    !!
    !! Creates a diagonal array or extract the diagonal elements of an array
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! diag-create-a-diagonal-array-or-extract-the-diagonal-elements-of-an-array))
      !
      ! Vector to matrix
      !
      procedure :: diag_rsp
      procedure :: diag_rdp
      procedure :: diag_csp
      procedure :: diag_cdp
      procedure :: diag_iint8
      procedure :: diag_iint16
      procedure :: diag_iint32
      procedure :: diag_iint64
      procedure :: diag_rsp_k
      procedure :: diag_rdp_k
      procedure :: diag_csp_k
      procedure :: diag_cdp_k
      procedure :: diag_iint8_k
      procedure :: diag_iint16_k
      procedure :: diag_iint32_k
      procedure :: diag_iint64_k

      !
      ! Matrix to vector
      !
      procedure :: diag_rsp_mat
      procedure :: diag_rdp_mat
      procedure :: diag_csp_mat
      procedure :: diag_cdp_mat
      procedure :: diag_iint8_mat
      procedure :: diag_iint16_mat
      procedure :: diag_iint32_mat
      procedure :: diag_iint64_mat
      procedure :: diag_rsp_mat_k
      procedure :: diag_rdp_mat_k
      procedure :: diag_csp_mat_k
      procedure :: diag_cdp_mat_k
      procedure :: diag_iint8_mat_k
      procedure :: diag_iint16_mat_k
      procedure :: diag_iint32_mat_k
      procedure :: diag_iint64_mat_k
  end interface


  ! Matrix trace
  interface trace
    !! version: experimental
    !!
    !! Computes the trace of a matrix
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! trace-trace-of-a-matrix))
      module procedure trace_rsp
      module procedure trace_rdp
      module procedure trace_csp
      module procedure trace_cdp
      module procedure trace_iint8
      module procedure trace_iint16
      module procedure trace_iint32
      module procedure trace_iint64
  end interface


  ! Outer product (of two vectors)
  interface outer_product
    !! version: experimental
    !!
    !! Computes the outer product of two vectors, returning a rank-2 array
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! outer_product-computes-the-outer-product-of-two-vectors))
      procedure :: outer_product_rsp
      procedure :: outer_product_rdp
      procedure :: outer_product_csp
      procedure :: outer_product_cdp
      procedure :: outer_product_iint8
      procedure :: outer_product_iint16
      procedure :: outer_product_iint32
      procedure :: outer_product_iint64
  end interface outer_product

  interface kronecker_product
    !! version: experimental
    !!
    !! Computes the Kronecker product of two arrays of size M1xN1, and of M2xN2, returning an (M1*M2)x(N1*N2) array
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! kronecker_product-computes-the-kronecker-product-of-two-matrices))
      procedure :: kronecker_product_rsp
      procedure :: kronecker_product_rdp
      procedure :: kronecker_product_csp
      procedure :: kronecker_product_cdp
      procedure :: kronecker_product_iint8
      procedure :: kronecker_product_iint16
      procedure :: kronecker_product_iint32
      procedure :: kronecker_product_iint64
  end interface kronecker_product


  ! Cross product (of two vectors)
  interface cross_product
    !! version: experimental
    !!
    !! Computes the cross product of two vectors, returning a rank-1 and size-3 array
    !! ([Specification](../page/specs/stdlib_linalg.html#cross_product-computes-the-cross-product-of-two-3-d-vectors))
      procedure :: cross_product_rsp
      procedure :: cross_product_rdp
      procedure :: cross_product_csp
      procedure :: cross_product_cdp
      procedure :: cross_product_iint8
      procedure :: cross_product_iint16
      procedure :: cross_product_iint32
      procedure :: cross_product_iint64
  end interface cross_product


  ! Check for squareness
  interface is_square
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is square
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_square-checks-if-a-matrix-is-square))
      module procedure is_square_rsp
      module procedure is_square_rdp
      module procedure is_square_csp
      module procedure is_square_cdp
      module procedure is_square_iint8
      module procedure is_square_iint16
      module procedure is_square_iint32
      module procedure is_square_iint64
  end interface is_square


  ! Check for diagonality
  interface is_diagonal
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is diagonal
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_diagonal-checks-if-a-matrix-is-diagonal))
      module procedure is_diagonal_rsp
      module procedure is_diagonal_rdp
      module procedure is_diagonal_csp
      module procedure is_diagonal_cdp
      module procedure is_diagonal_iint8
      module procedure is_diagonal_iint16
      module procedure is_diagonal_iint32
      module procedure is_diagonal_iint64
  end interface is_diagonal


  ! Check for symmetry
  interface is_symmetric
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is symmetric
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_symmetric-checks-if-a-matrix-is-symmetric))
      module procedure is_symmetric_rsp
      module procedure is_symmetric_rdp
      module procedure is_symmetric_csp
      module procedure is_symmetric_cdp
      module procedure is_symmetric_iint8
      module procedure is_symmetric_iint16
      module procedure is_symmetric_iint32
      module procedure is_symmetric_iint64
  end interface is_symmetric


  ! Check for skew-symmetry
  interface is_skew_symmetric
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is skew-symmetric
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_skew_symmetric-checks-if-a-matrix-is-skew-symmetric))
      module procedure is_skew_symmetric_rsp
      module procedure is_skew_symmetric_rdp
      module procedure is_skew_symmetric_csp
      module procedure is_skew_symmetric_cdp
      module procedure is_skew_symmetric_iint8
      module procedure is_skew_symmetric_iint16
      module procedure is_skew_symmetric_iint32
      module procedure is_skew_symmetric_iint64
  end interface is_skew_symmetric


  ! Check for Hermiticity
  interface is_hermitian
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is Hermitian
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_hermitian-checks-if-a-matrix-is-hermitian))
      module procedure is_hermitian_rsp
      module procedure is_hermitian_rdp
      module procedure is_hermitian_csp
      module procedure is_hermitian_cdp
      module procedure is_hermitian_iint8
      module procedure is_hermitian_iint16
      module procedure is_hermitian_iint32
      module procedure is_hermitian_iint64
  end interface is_hermitian


  ! Check for triangularity
  interface is_triangular
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is triangular
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_triangular-checks-if-a-matrix-is-triangular))
      module procedure is_triangular_rsp
      module procedure is_triangular_rdp
      module procedure is_triangular_csp
      module procedure is_triangular_cdp
      module procedure is_triangular_iint8
      module procedure is_triangular_iint16
      module procedure is_triangular_iint32
      module procedure is_triangular_iint64
  end interface is_triangular
  

  ! Check for matrix being Hessenberg
  interface is_hessenberg
    !! version: experimental
    !!
    !! Checks if a matrix (rank-2 array) is Hessenberg
    !! ([Specification](../page/specs/stdlib_linalg.html#
    !! is_hessenberg-checks-if-a-matrix-is-hessenberg))
      module procedure is_Hessenberg_rsp
      module procedure is_Hessenberg_rdp
      module procedure is_Hessenberg_csp
      module procedure is_Hessenberg_cdp
      module procedure is_Hessenberg_iint8
      module procedure is_Hessenberg_iint16
      module procedure is_Hessenberg_iint32
      module procedure is_Hessenberg_iint64
  end interface is_hessenberg

contains


    !> Version: experimental
    !>
    !> Constructs the identity matrix.
    !> ([Specification](../page/specs/stdlib_linalg.html#eye-construct-the-identity-matrix))
    pure function eye(dim1, dim2) result(result)

        integer, intent(in) :: dim1
        integer, intent(in), optional :: dim2
        integer(int8), allocatable :: result(:, :)

        integer :: dim2_
        integer :: i

        dim2_ = optval(dim2, dim1)
        allocate(result(dim1, dim2_))
        
        result = 0_int8
        do i = 1, min(dim1, dim2_)
            result(i, i) = 1_int8
        end do

    end function eye

      function trace_rsp(A) result(res)
        real(sp), intent(in) :: A(:,:)
        real(sp) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_rsp
      function trace_rdp(A) result(res)
        real(dp), intent(in) :: A(:,:)
        real(dp) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_rdp
      function trace_csp(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        complex(sp) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_csp
      function trace_cdp(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        complex(dp) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_cdp
      function trace_iint8(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        integer(int8) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_iint8
      function trace_iint16(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        integer(int16) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_iint16
      function trace_iint32(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        integer(int32) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_iint32
      function trace_iint64(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        integer(int64) :: res
        integer :: i
        res = 0
        do i = 1, minval(shape(A))
          res = res + A(i,i)
        end do
      end function trace_iint64


      pure function is_square_rsp(A) result(res)
        real(sp), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_rsp
      pure function is_square_rdp(A) result(res)
        real(dp), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_rdp
      pure function is_square_csp(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_csp
      pure function is_square_cdp(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_cdp
      pure function is_square_iint8(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_iint8
      pure function is_square_iint16(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_iint16
      pure function is_square_iint32(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_iint32
      pure function is_square_iint64(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        logical :: res
        res = (size(A,1) == size(A,2))
      end function is_square_iint64


      pure function is_diagonal_rsp(A) result(res)
        real(sp), intent(in) :: A(:,:)
        logical :: res
        real(sp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_rsp
      pure function is_diagonal_rdp(A) result(res)
        real(dp), intent(in) :: A(:,:)
        logical :: res
        real(dp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_rdp
      pure function is_diagonal_csp(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        logical :: res
        complex(sp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_csp
      pure function is_diagonal_cdp(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        logical :: res
        complex(dp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_cdp
      pure function is_diagonal_iint8(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        logical :: res
        integer(int8), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_iint8
      pure function is_diagonal_iint16(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        logical :: res
        integer(int16), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_iint16
      pure function is_diagonal_iint32(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        logical :: res
        integer(int32), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_iint32
      pure function is_diagonal_iint64(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        logical :: res
        integer(int64), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        do j = 1, n !loop over all columns
            o = min(j-1,m) !index of row above diagonal (or last row)
            do i = 1, o !loop over rows above diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
            do i = o+2, m !loop over rows below diagonal
                if (A(i,j) /= zero) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is diagonal
      end function is_diagonal_iint64


      pure function is_symmetric_rsp(A) result(res)
        real(sp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_rsp
      pure function is_symmetric_rdp(A) result(res)
        real(dp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_rdp
      pure function is_symmetric_csp(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_csp
      pure function is_symmetric_cdp(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_cdp
      pure function is_symmetric_iint8(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_iint8
      pure function is_symmetric_iint16(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_iint16
      pure function is_symmetric_iint32(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_iint32
      pure function is_symmetric_iint64(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j-1 !loop over all rows above diagonal
                if (A(i,j) /= A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is symmetric
      end function is_symmetric_iint64


      pure function is_skew_symmetric_rsp(A) result(res)
        real(sp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_rsp
      pure function is_skew_symmetric_rdp(A) result(res)
        real(dp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_rdp
      pure function is_skew_symmetric_csp(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_csp
      pure function is_skew_symmetric_cdp(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_cdp
      pure function is_skew_symmetric_iint8(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_iint8
      pure function is_skew_symmetric_iint16(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_iint16
      pure function is_skew_symmetric_iint32(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_iint32
      pure function is_skew_symmetric_iint64(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be skew-symmetric
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= -A(j,i)) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is skew-symmetric
      end function is_skew_symmetric_iint64


      pure function is_hermitian_rsp(A) result(res)
        real(sp), intent(in) :: A(:,:)
        logical :: res
        res = is_symmetric(A) !symmetry and Hermiticity are equivalent for real matrices
      end function is_hermitian_rsp
      pure function is_hermitian_rdp(A) result(res)
        real(dp), intent(in) :: A(:,:)
        logical :: res
        res = is_symmetric(A) !symmetry and Hermiticity are equivalent for real matrices
      end function is_hermitian_rdp
      pure function is_hermitian_iint8(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        logical :: res
        res = is_symmetric(A) !symmetry and Hermiticity are equivalent for real matrices
      end function is_hermitian_iint8
      pure function is_hermitian_iint16(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        logical :: res
        res = is_symmetric(A) !symmetry and Hermiticity are equivalent for real matrices
      end function is_hermitian_iint16
      pure function is_hermitian_iint32(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        logical :: res
        res = is_symmetric(A) !symmetry and Hermiticity are equivalent for real matrices
      end function is_hermitian_iint32
      pure function is_hermitian_iint64(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        logical :: res
        res = is_symmetric(A) !symmetry and Hermiticity are equivalent for real matrices
      end function is_hermitian_iint64
      pure function is_hermitian_csp(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be Hermitian
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= conjg(A(j,i))) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is Hermitian
      end function is_hermitian_csp
      pure function is_hermitian_cdp(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        logical :: res
        integer :: n, i, j
        if (.not. is_square(A)) then
           res = .false.
           return !nonsquare matrices cannot be Hermitian
        end if
        n = size(A,1) !symmetric dimension of A
        do j = 1, n !loop over all columns
            do i = 1, j !loop over all rows above diagonal (and diagonal)
                if (A(i,j) /= conjg(A(j,i))) then
                  res = .false.
                  return
                end if
            end do
        end do
        res = .true. !otherwise A is Hermitian
      end function is_hermitian_cdp


      function is_triangular_rsp(A,uplo) result(res)
        real(sp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        real(sp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_rsp
      function is_triangular_rdp(A,uplo) result(res)
        real(dp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        real(dp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_rdp
      function is_triangular_csp(A,uplo) result(res)
        complex(sp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        complex(sp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_csp
      function is_triangular_cdp(A,uplo) result(res)
        complex(dp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        complex(dp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_cdp
      function is_triangular_iint8(A,uplo) result(res)
        integer(int8), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int8), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_iint8
      function is_triangular_iint16(A,uplo) result(res)
        integer(int16), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int16), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_iint16
      function is_triangular_iint32(A,uplo) result(res)
        integer(int32), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int32), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_iint32
      function is_triangular_iint64(A,uplo) result(res)
        integer(int64), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int64), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper triangularity
          do j = 1, n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i = o+2, m !loop over rows below diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower triangularity
          do j=1,n !loop over all columns
              o = min(j-1,m) !index of row above diagonal (or last row)
              do i=1,o !loop over rows above diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
     
        res = .true. !otherwise A is triangular of the requested type
      end function is_triangular_iint64


      function is_hessenberg_rsp(A,uplo) result(res)
        real(sp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        real(sp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_rsp
      function is_hessenberg_rdp(A,uplo) result(res)
        real(dp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        real(dp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_rdp
      function is_hessenberg_csp(A,uplo) result(res)
        complex(sp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        complex(sp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_csp
      function is_hessenberg_cdp(A,uplo) result(res)
        complex(dp), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        complex(dp), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_cdp
      function is_hessenberg_iint8(A,uplo) result(res)
        integer(int8), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int8), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_iint8
      function is_hessenberg_iint16(A,uplo) result(res)
        integer(int16), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int16), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_iint16
      function is_hessenberg_iint32(A,uplo) result(res)
        integer(int32), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int32), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_iint32
      function is_hessenberg_iint64(A,uplo) result(res)
        integer(int64), intent(in) :: A(:,:)
        character, intent(in) :: uplo
        logical :: res
        integer(int64), parameter :: zero = 0 !zero of relevant type
        integer :: m, n, o, i, j
        m = size(A,1)
        n = size(A,2)
        if ((uplo == 'u') .or. (uplo == 'U')) then !check for upper Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = o+4, m !loop over rows two or more below main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
          end do
        else if ((uplo == 'l') .or. (uplo == 'L')) then !check for lower Hessenberg
          do j = 1, n !loop over all columns
              o = min(j-2,m) !index of row two above diagonal (or last row)
              do i = 1, o !loop over rows one or more above main diagonal
                  if (A(i,j) /= zero) then
                    res = .false.
                    return
                  end if
              end do
           end do
        else
           error stop
        end if
        res = .true. !otherwise A is Hessenberg of the requested type
      end function is_hessenberg_iint64
    

      pure function diag_rsp(v) result(res)
        real(sp), intent(in) :: v(:)
        real(sp) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_rsp
      pure function diag_rdp(v) result(res)
        real(dp), intent(in) :: v(:)
        real(dp) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_rdp
      pure function diag_csp(v) result(res)
        complex(sp), intent(in) :: v(:)
        complex(sp) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_csp
      pure function diag_cdp(v) result(res)
        complex(dp), intent(in) :: v(:)
        complex(dp) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_cdp
      pure function diag_iint8(v) result(res)
        integer(int8), intent(in) :: v(:)
        integer(int8) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_iint8
      pure function diag_iint16(v) result(res)
        integer(int16), intent(in) :: v(:)
        integer(int16) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_iint16
      pure function diag_iint32(v) result(res)
        integer(int32), intent(in) :: v(:)
        integer(int32) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_iint32
      pure function diag_iint64(v) result(res)
        integer(int64), intent(in) :: v(:)
        integer(int64) :: res(size(v),size(v))
        integer :: i
        res = 0
        do i = 1, size(v)
          res(i,i) = v(i)
        end do
      end function diag_iint64


      pure function diag_rsp_k(v,k) result(res)
        real(sp), intent(in) :: v(:)
        integer, intent(in) :: k
        real(sp) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_rsp_k
      pure function diag_rdp_k(v,k) result(res)
        real(dp), intent(in) :: v(:)
        integer, intent(in) :: k
        real(dp) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_rdp_k
      pure function diag_csp_k(v,k) result(res)
        complex(sp), intent(in) :: v(:)
        integer, intent(in) :: k
        complex(sp) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_csp_k
      pure function diag_cdp_k(v,k) result(res)
        complex(dp), intent(in) :: v(:)
        integer, intent(in) :: k
        complex(dp) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_cdp_k
      pure function diag_iint8_k(v,k) result(res)
        integer(int8), intent(in) :: v(:)
        integer, intent(in) :: k
        integer(int8) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_iint8_k
      pure function diag_iint16_k(v,k) result(res)
        integer(int16), intent(in) :: v(:)
        integer, intent(in) :: k
        integer(int16) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_iint16_k
      pure function diag_iint32_k(v,k) result(res)
        integer(int32), intent(in) :: v(:)
        integer, intent(in) :: k
        integer(int32) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_iint32_k
      pure function diag_iint64_k(v,k) result(res)
        integer(int64), intent(in) :: v(:)
        integer, intent(in) :: k
        integer(int64) :: res(size(v)+abs(k),size(v)+abs(k))
        integer :: i, sz
        sz = size(v)
        res = 0
        if (k > 0) then
          do i = 1, sz
              res(i,k+i) = v(i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i+abs(k),i) = v(i)
          end do
        else
          do i = 1, sz
              res(i,i) = v(i)
          end do
        end if
      end function diag_iint64_k

      pure function diag_rsp_mat(A) result(res)
        real(sp), intent(in) :: A(:,:)
        real(sp) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_rsp_mat
      pure function diag_rdp_mat(A) result(res)
        real(dp), intent(in) :: A(:,:)
        real(dp) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_rdp_mat
      pure function diag_csp_mat(A) result(res)
        complex(sp), intent(in) :: A(:,:)
        complex(sp) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_csp_mat
      pure function diag_cdp_mat(A) result(res)
        complex(dp), intent(in) :: A(:,:)
        complex(dp) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_cdp_mat
      pure function diag_iint8_mat(A) result(res)
        integer(int8), intent(in) :: A(:,:)
        integer(int8) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_iint8_mat
      pure function diag_iint16_mat(A) result(res)
        integer(int16), intent(in) :: A(:,:)
        integer(int16) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_iint16_mat
      pure function diag_iint32_mat(A) result(res)
        integer(int32), intent(in) :: A(:,:)
        integer(int32) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_iint32_mat
      pure function diag_iint64_mat(A) result(res)
        integer(int64), intent(in) :: A(:,:)
        integer(int64) :: res(minval(shape(A)))
        integer :: i
        do i = 1, minval(shape(A))
          res(i) = A(i,i)
        end do
      end function diag_iint64_mat

      pure function diag_rsp_mat_k(A,k) result(res)
        real(sp), intent(in) :: A(:,:)
        integer, intent(in) :: k
        real(sp) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_rsp_mat_k
      pure function diag_rdp_mat_k(A,k) result(res)
        real(dp), intent(in) :: A(:,:)
        integer, intent(in) :: k
        real(dp) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_rdp_mat_k
      pure function diag_csp_mat_k(A,k) result(res)
        complex(sp), intent(in) :: A(:,:)
        integer, intent(in) :: k
        complex(sp) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_csp_mat_k
      pure function diag_cdp_mat_k(A,k) result(res)
        complex(dp), intent(in) :: A(:,:)
        integer, intent(in) :: k
        complex(dp) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_cdp_mat_k
      pure function diag_iint8_mat_k(A,k) result(res)
        integer(int8), intent(in) :: A(:,:)
        integer, intent(in) :: k
        integer(int8) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_iint8_mat_k
      pure function diag_iint16_mat_k(A,k) result(res)
        integer(int16), intent(in) :: A(:,:)
        integer, intent(in) :: k
        integer(int16) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_iint16_mat_k
      pure function diag_iint32_mat_k(A,k) result(res)
        integer(int32), intent(in) :: A(:,:)
        integer, intent(in) :: k
        integer(int32) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_iint32_mat_k
      pure function diag_iint64_mat_k(A,k) result(res)
        integer(int64), intent(in) :: A(:,:)
        integer, intent(in) :: k
        integer(int64) :: res(minval(shape(A))-abs(k))
        integer :: i, sz
        sz = minval(shape(A))-abs(k)
        if (k > 0) then
          do i = 1, sz
              res(i) = A(i,k+i)
          end do
        else if (k < 0) then
          do i = 1, sz
              res(i) = A(i+abs(k),i)
          end do
        else
          do i = 1, sz
              res(i) = A(i,i)
          end do
        end if
      end function diag_iint64_mat_k

    pure function outer_product_rsp(u, v) result(res)
      real(sp), intent(in) :: u(:), v(:)
      real(sp) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_rsp
    pure function outer_product_rdp(u, v) result(res)
      real(dp), intent(in) :: u(:), v(:)
      real(dp) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_rdp
    pure function outer_product_csp(u, v) result(res)
      complex(sp), intent(in) :: u(:), v(:)
      complex(sp) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_csp
    pure function outer_product_cdp(u, v) result(res)
      complex(dp), intent(in) :: u(:), v(:)
      complex(dp) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_cdp
    pure function outer_product_iint8(u, v) result(res)
      integer(int8), intent(in) :: u(:), v(:)
      integer(int8) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_iint8
    pure function outer_product_iint16(u, v) result(res)
      integer(int16), intent(in) :: u(:), v(:)
      integer(int16) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_iint16
    pure function outer_product_iint32(u, v) result(res)
      integer(int32), intent(in) :: u(:), v(:)
      integer(int32) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_iint32
    pure function outer_product_iint64(u, v) result(res)
      integer(int64), intent(in) :: u(:), v(:)
      integer(int64) :: res(size(u),size(v))
      integer :: col
      do col = 1, size(v)
        res(:,col) = v(col) * u
      end do
    end function outer_product_iint64


    pure function kronecker_product_rsp(A, B) result(C)
      real(sp), intent(in) :: A(:,:), B(:,:)
      real(sp) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_rsp
    pure function kronecker_product_rdp(A, B) result(C)
      real(dp), intent(in) :: A(:,:), B(:,:)
      real(dp) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_rdp
    pure function kronecker_product_csp(A, B) result(C)
      complex(sp), intent(in) :: A(:,:), B(:,:)
      complex(sp) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_csp
    pure function kronecker_product_cdp(A, B) result(C)
      complex(dp), intent(in) :: A(:,:), B(:,:)
      complex(dp) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_cdp
    pure function kronecker_product_iint8(A, B) result(C)
      integer(int8), intent(in) :: A(:,:), B(:,:)
      integer(int8) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_iint8
    pure function kronecker_product_iint16(A, B) result(C)
      integer(int16), intent(in) :: A(:,:), B(:,:)
      integer(int16) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_iint16
    pure function kronecker_product_iint32(A, B) result(C)
      integer(int32), intent(in) :: A(:,:), B(:,:)
      integer(int32) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_iint32
    pure function kronecker_product_iint64(A, B) result(C)
      integer(int64), intent(in) :: A(:,:), B(:,:)
      integer(int64) :: C(size(A,dim=1)*size(B,dim=1),size(A,dim=2)*size(B,dim=2))
      integer :: m1, n1, maxM1, maxN1, maxM2, maxN2
      
      maxM1 = size(A, dim=1)
      maxN1 = size(A, dim=2)
      maxM2 = size(B, dim=1)
      maxN2 = size(B, dim=2)
      

      do n1 = 1, maxN1
         do m1 = 1, maxM1
            ! We use the Wikipedia convention for ordering of the matrix elements
	    ! https://en.wikipedia.org/wiki/Kronecker_product
            C((m1-1)*maxM2+1:m1*maxM2, (n1-1)*maxN2+1:n1*maxN2) = A(m1, n1) * B(:,:)
         end do
      end do
    end function kronecker_product_iint64

  pure function cross_product_rsp(a, b) result(res)
      real(sp), intent(in) :: a(3), b(3)
      real(sp) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_rsp
  pure function cross_product_rdp(a, b) result(res)
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_rdp
  pure function cross_product_csp(a, b) result(res)
      complex(sp), intent(in) :: a(3), b(3)
      complex(sp) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_csp
  pure function cross_product_cdp(a, b) result(res)
      complex(dp), intent(in) :: a(3), b(3)
      complex(dp) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_cdp
  pure function cross_product_iint8(a, b) result(res)
      integer(int8), intent(in) :: a(3), b(3)
      integer(int8) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_iint8
  pure function cross_product_iint16(a, b) result(res)
      integer(int16), intent(in) :: a(3), b(3)
      integer(int16) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_iint16
  pure function cross_product_iint32(a, b) result(res)
      integer(int32), intent(in) :: a(3), b(3)
      integer(int32) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_iint32
  pure function cross_product_iint64(a, b) result(res)
      integer(int64), intent(in) :: a(3), b(3)
      integer(int64) :: res(3)

      res(1) = a(2) * b(3) - a(3) * b(2)
      res(2) = a(3) * b(1) - a(1) * b(3)
      res(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product_iint64
end module stdlib_linalg
