module mbd_util

use mbd_interface, only: print_error, print_log

implicit none

private

public :: tostr, start_clock, stop_clock, get_plock, release_plock, &
    have_plock, operator(.cprod.), diag, eye, invert, inverted, &
    diagonalize, sdiagonalize, diagonalized, sdiagonalized, solve_linear_eqs

interface tostr
    module procedure tostr_int_
    module procedure tostr_dble_
end interface

type :: Timestamp_t
    character(len=30) :: label = ''
    integer :: cnt = 0  ! number of cpu clocks
    integer :: n = 0  ! how many times was run
end type

character(len=50) :: parallel_lock = ''

logical, public :: measure_time = .false.
integer, parameter :: n_timestamps = 50
type(Timestamp_t), save :: timestamps(n_timestamps)

interface operator(.cprod.)
    module procedure cart_prod_
end interface

interface diag
    module procedure get_diag_
    module procedure get_diag_cmplx_
    module procedure make_diag_
end interface

interface invert
    module procedure invert_ge_dble_
    module procedure invert_ge_cmplx_
end interface

interface diagonalize
    module procedure diagonalize_ge_dble_
    module procedure diagonalize_ge_cmplx_
end interface

interface sdiagonalize
    module procedure diagonalize_sym_dble_
    module procedure diagonalize_he_cmplx_
end interface

interface diagonalized
    module procedure diagonalized_ge_dble_
end interface

interface sdiagonalized
    module procedure diagonalized_sym_dble_
end interface

external :: ZHEEV, DGEEV, DSYEV, DGETRF, DGETRI, DGESV, ZGETRF, ZGETRI, &
    ZGEEV, ZGEEB

contains

integer function get_stamp(label) result(i_stamp)
    character(len=*), intent(in) :: label

    do i_stamp = 1, n_timestamps
        if (timestamps(i_stamp)%label == label) then
            exit
        else if (timestamps(i_stamp)%label == '') then
            timestamps(i_stamp)%label = label
            exit
        end if
    end do
end function

subroutine start_clock(label)
    character(len=*), intent(in) :: label

    integer :: cnt, rate, cnt_max, i_stamp

    if (.not. measure_time) return
    i_stamp = get_stamp(label)
    if (i_stamp > n_timestamps) return
    call system_clock(cnt, rate, cnt_max)
    timestamps(i_stamp)%cnt = timestamps(i_stamp)%cnt - cnt
    timestamps(i_stamp)%n= timestamps(i_stamp)%n + 1
end subroutine

subroutine stop_clock(label)
    character(len=*), intent(in) :: label

    integer :: cnt, rate, cnt_max, i_stamp

    if (.not. measure_time) return
    i_stamp = get_stamp(label)
    if (i_stamp > n_timestamps) return
    call system_clock(cnt, rate, cnt_max)
    timestamps(i_stamp)%cnt = timestamps(i_Stamp)%cnt + cnt
end subroutine

subroutine get_plock(label)
    character(len=*), intent(in) :: label

    if (parallel_lock == '') then
        parallel_lock = label
    end if
end subroutine

logical function have_plock(label)
    character(len=*), intent(in) :: label

    have_plock = parallel_lock == label
end function

subroutine release_plock(label)
    character(len=*), intent(in) :: label

    if (parallel_lock == label) then
        parallel_lock = ''
    end if
end subroutine

character(len=50) elemental function tostr_int_(k, format)
    implicit none

    integer, intent(in) :: k
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_int_, format) k
    else
        write (tostr_int_, "(i20)") k
    end if
    tostr_int_ = adjustl(tostr_int_)
end function tostr_int_

character(len=50) elemental function tostr_dble_(x, format)
    implicit none

    double precision, intent(in) :: x
    character(*), intent(in), optional :: format

    if (present(format)) then
        write (tostr_dble_, format) x
    else
        write (tostr_dble_, "(g50.17e3)") x
    end if
    tostr_dble_ = adjustl(tostr_dble_)
end function tostr_dble_

function eye(n) result(A)
    integer, intent(in) :: n
    real(8) :: A(n, n)

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:n) A(i, i) = 1.d0
end function

subroutine invert_ge_dble_(A)
    real(8), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    real(8), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    real(8) :: n_work_arr_optim
    integer :: error_flag

    call start_clock('inversion')
    n = size(A, 1)
    call DGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            'Matrix inversion failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    call DGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim)
    allocate (work_arr(n_work_arr))
    call DGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_error( &
            'Matrix inversion failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    call stop_clock('inversion')
end subroutine

subroutine invert_ge_cmplx_(A)
    complex(8), intent(inout) :: A(:, :)

    integer :: i_pivot(size(A, 1))
    complex(8), allocatable :: work_arr(:)
    integer :: n
    integer :: n_work_arr
    complex(8) :: n_work_arr_optim
    integer :: error_flag

    call start_clock('inversion')
    n = size(A, 1)
    call ZGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        call print_error( &
            'Matrix inversion failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    call ZGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(dble(n_work_arr_optim))
    allocate (work_arr(n_work_arr))
    call ZGETRI(n, A, n, i_pivot, work_arr, n_work_arr, error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_error( &
            'Matrix inversion failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    call stop_clock('inversion')
end subroutine

function inverted(A) result(A_inv)
    real(8), intent(in) :: A(:, :)
    real(8) :: A_inv(size(A, 1), size(A, 2))

    A_inv = A
    call invert(A_inv)
end function

subroutine diagonalize_sym_dble_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    real(8), intent(inout) :: A(:, :)
    real(8), intent(out) :: eigs(:)

    real(8), allocatable :: work_arr(:)
    integer :: n
    real(8) :: n_work_arr
    integer :: error_flag

    call start_clock('diagonalization')
    n = size(A, 1)
    call DSYEV(mode, 'U', n, A, n, eigs, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DSYEV(mode, 'U', n, A, n, eigs, work_arr, size(work_arr), error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_log( &
            'DSYEV failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    call stop_clock('diagonalization')
end subroutine

function diagonalized_sym_dble_(A) result(eigs)
    real(8), intent(in) :: A(:, :)
    real(8) :: eigs(size(A, 1))

    real(8) :: eigvecs(size(A, 1), size(A, 2))

    eigvecs = A
    call sdiagonalize('N', eigvecs, eigs)
end function

subroutine diagonalize_ge_dble_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    real(8), intent(inout) :: A(:, :)
    complex(8), intent(out) :: eigs(size(A, 1))

    real(8), allocatable :: work_arr(:)
    integer :: n
    real(8) :: n_work_arr
    integer :: error_flag
    real(8) :: eigs_r(size(A, 1)), eigs_i(size(A, 1))
    real(8) :: dummy
    real(8) :: vectors(size(A, 1), size(A, 2))

    call start_clock('diagonalization')
    n = size(A, 1)
    call DGEEV('N', mode, n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr)))
    call DGEEV('N', mode, n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, work_arr, size(work_arr), error_flag)
    deallocate (work_arr)
    if (error_flag /= 0) then
        call print_log( &
            'DGEEV failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    eigs = cmplx(eigs_r, eigs_i, 8)
    A = vectors
    call stop_clock('diagonalization')
end subroutine

function diagonalized_ge_dble_(A) result(eigs)
    real(8), intent(in) :: A(:, :)
    complex(8) :: eigs(size(A, 1))

    real(8) :: eigvecs(size(A, 1), size(A, 2))

    eigvecs = A
    call diagonalize('N', eigvecs, eigs)
end function

subroutine diagonalize_he_cmplx_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    complex(8), intent(inout) :: A(:, :)
    real(8), intent(out) :: eigs(size(A, 1))

    complex(8), allocatable :: work(:)
    complex(8) :: lwork_cmplx
    real(8), allocatable :: rwork(:)
    integer :: n, lwork
    integer :: error_flag
    integer, external :: ILAENV

    call start_clock('diagonalization')
    n = size(A, 1)
    allocate (rwork(max(1, 3*n-2)))
    call ZHEEV(mode, 'U', n, A, n, eigs, lwork_cmplx, -1, rwork, error_flag)
    lwork = nint(dble(lwork_cmplx))
    allocate (work(lwork))
    call ZHEEV(mode, 'U', n, A, n, eigs, work, lwork, rwork, error_flag)
    deallocate (rwork)
    deallocate (work)
    if (error_flag /= 0) then
        call print_error( &
            'ZHEEV failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    call stop_clock('diagonalization')
end subroutine

subroutine diagonalize_ge_cmplx_(mode, A, eigs)
    character(len=1), intent(in) :: mode
    complex(8), intent(inout) :: A(:, :)
    complex(8), intent(out) :: eigs(size(A, 1))

    complex(8), allocatable :: work(:)
    real(8) :: rwork(2*size(A, 1))
    integer :: n, lwork
    complex(8) :: lwork_arr
    integer :: error_flag
    complex(8) :: dummy
    complex(8) :: vectors(size(A, 1), size(A, 2))

    call start_clock('diagonalization')
    n = size(A, 1)
    call ZGEEV('N', mode, n, A, n, eigs, dummy, 1, &
        vectors, n, lwork_arr, -1, rwork, error_flag)
    lwork = nint(dble(lwork_arr))
    allocate (work(lwork))
    call ZGEEV('N', mode, n, A, n, eigs, dummy, 1, &
        vectors, n, work, lwork, rwork, error_flag)
    deallocate (work)
    if (error_flag /= 0) then
        call print_log( &
            'ZGEEV failed in module mbd with error code ' &
            // trim(tostr(error_flag)) &
        )
    endif
    A = vectors
    call stop_clock('diagonalization')
end subroutine

function cart_prod_(a, b) result(c)
    real(8), intent(in) :: a(:), b(:)
    real(8) :: c(size(a), size(b))

    integer :: i, j

    forall (i = 1:size(a), j = 1:size(b)) c(i, j) = a(i)*b(j)
end function

function get_diag_(A) result(d)
    real(8), intent(in) :: A(:, :)
    real(8) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function

function get_diag_cmplx_(A) result(d)
    complex(8), intent(in) :: A(:, :)
    complex(8) :: d(size(A, 1))

    integer :: i

    forall (i = 1:size(A, 1)) d(i) = A(i, i)
end function

function make_diag_(d) result(A)
    real(8), intent(in) :: d(:)
    real(8) :: A(size(d), size(d))

    integer :: i

    A(:, :) = 0.d0
    forall (i = 1:size(d)) A(i, i) = d(i)
end function

function solve_linear_eqs(A, b) result(x)
    real(8), intent(in) :: A(:, :), b(:)
    real(8) :: x(size(b))

    real(8) :: A_(size(b), size(b))
    integer :: i_pivot(size(b))
    integer :: n
    integer :: error_flag

    A_ = A
    x = b
    n = size(b)
    call DGESV(n, 1, A_, n, i_pivot, x, n, error_flag)
end function

end module
