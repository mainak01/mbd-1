module mbd_interface

use iso_c_binding
use mpi

implicit none

private

integer(c_int), public, bind(c) :: my_task = 0, n_tasks = 1

public :: &
    sync_sum, print_log, print_error, print_warning, &
    mute, unmute

integer :: muted = 0

interface sync_sum
    module procedure sync_sum_dble_
    module procedure sync_sum_vector_dble_
    module procedure sync_sum_matrix_dble_
    module procedure sync_sum_3d_dble_
    module procedure sync_sum_4d_dble_
    module procedure sync_sum_cplx_
    module procedure sync_sum_vector_cplx_
    module procedure sync_sum_matrix_cplx_
    module procedure sync_sum_3d_cplx_
    module procedure sync_sum_4d_cplx_
end interface

external :: MPI_Allreduce

contains

subroutine sync_sum_dble_(x)
    real(8), intent(inout) :: x
    real(8) :: x_buff
    integer :: mpi_err

    call MPI_Allreduce( &
        x, x_buff, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    x = x_buff
end subroutine

subroutine sync_sum_array_dble_(array, n_array)
    integer, intent(in) :: n_array
    real(8), intent(inout) :: array(n_array)

    real(8) :: array_buff(n_array)
    integer :: mpi_err

    call MPI_Allreduce( &
        array, array_buff, n_array, MPI_DOUBLE_PRECISION, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    array = array_buff
end subroutine

subroutine sync_sum_vector_dble_(x)
    real(8), intent(inout) :: x(:)

    call sync_sum_array_dble_(x, size(x))
end subroutine

subroutine sync_sum_matrix_dble_(x)
    real(8), intent(inout) :: x(:, :)

    call  sync_sum_array_dble_(x, size(x))
end subroutine

subroutine sync_sum_3d_dble_(x)
    real(8), intent(inout) :: x(:, :, :)

    call  sync_sum_array_dble_(x, size(x))
end subroutine

subroutine sync_sum_4d_dble_(x)
    real(8), intent(inout) :: x(:, :, :, :)

    call  sync_sum_array_dble_(x, size(x))
end subroutine

subroutine sync_sum_cplx_(x)
    complex(kind=8), intent(inout) :: x
    complex(kind=8) :: x_buff
    integer :: mpi_err

    call MPI_Allreduce( &
        x, x_buff, 1, MPI_COMPLEX16, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    x = x_buff
end subroutine

subroutine sync_sum_array_cplx_(array, n_array)
    integer, intent(in) :: n_array
    complex(kind=8), intent(inout) :: array(n_array)

    complex(kind=8) :: array_buff(n_array)
    integer :: mpi_err

    call MPI_Allreduce( &
        array, array_buff, n_array, MPI_COMPLEX16, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)
    array = array_buff
end subroutine

subroutine sync_sum_vector_cplx_(x)
    complex(kind=8), intent(inout) :: x(:)

    call sync_sum_array_cplx_(x, size(x))
end subroutine

subroutine sync_sum_matrix_cplx_(x)
    complex(kind=8), intent(inout) :: x(:, :)

    call  sync_sum_array_cplx_(x, size(x))
end subroutine

subroutine sync_sum_3d_cplx_(x)
    complex(kind=8), intent(inout) :: x(:, :, :)

    call  sync_sum_array_cplx_(x, size(x))
end subroutine

subroutine sync_sum_4d_cplx_(x)
    complex(kind=8), intent(inout) :: x(:, :, :, :)

    call  sync_sum_array_cplx_(x, size(x))
end subroutine

subroutine mute()
    muted = muted + 1
end subroutine

subroutine unmute()
    muted = muted - 1
end subroutine

subroutine print_log(str)
    character(len=*), intent(in) :: str

    integer :: myid, error

    if (muted > 0) return
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (6, *) str
    end if
end subroutine

subroutine print_warning(str)
    character(len=*), intent(in) :: str

    integer :: myid, error

    call MPI_Comm_rank(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (0, *) "Warning: " // str
    end if
end subroutine

subroutine print_error(str)
    character(len=*), intent(in) :: str

    integer :: myid, error

    call MPI_Comm_rank(MPI_COMM_WORLD, myid, error)
    if (myid == 0) then
        write (0, *) "Error: " // str
    end if
end subroutine

end module mbd_interface
