!MPI helpers for the distributed unit tests (compiled only with mpi=yes).
!
!The process grid is given on the command line as `npy npz`, it is the 2D
!PFFT process mesh: x is never distributed, y is split into npy and z
!into npz parts. The local block of the global grid is obtained from
!PoisFFT_LocalGridSize, i.e. the same decomposition the solvers expect.
module mpi_test_grid
  use mpi
  use iso_c_binding, only: c_int, c_intptr_t
  use iso_fortran_env, only: dp => real64, error_unit
  use ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use PoisFFT, only: PoisFFT_InitMPIGrid, PoisFFT_LocalGridSize
  implicit none

  private
  public :: grid_init, grid_finish, global_rel_err, cart_comm, rank, nprocs, &
            npy, npz, nxyz, offs, is_root, command_case

  integer :: cart_comm = MPI_COMM_NULL
  integer :: rank = 0, nprocs = 1
  integer :: npy = 1, npz = 1
  integer :: nxyz(3), offs(3)
  logical :: is_root = .true.
  !the third command line argument, selects the test case
  integer :: command_case = 0

contains

  subroutine grid_init(gnxyz)
    integer, intent(in) :: gnxyz(3)
    integer :: ie, provided
    character(32) :: arg
    integer(c_intptr_t) :: ln(3), lo(3), ln2(3), lo2(3)

    call MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ie)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ie)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ie)
    is_root = rank == 0

    !the solvers call MPI from inside OpenMP single regions
    !$ if (provided < MPI_THREAD_SERIALIZED) then
    !$   if (is_root) write(error_unit, '(a)') "MPI does not provide MPI_THREAD_SERIALIZED."
    !$   call MPI_Abort(MPI_COMM_WORLD, 5, ie)
    !$ end if

    if (command_argument_count() >= 2) then
      call get_command_argument(1, arg)
      read(arg, *) npy
      call get_command_argument(2, arg)
      read(arg, *) npz
    else
      npy = 1
      npz = nprocs
    end if
    if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read(arg, *) command_case
    end if

    if (npy * npz /= nprocs) then
      if (is_root) write(error_unit, '(a,i0,a,i0,a,i0)') &
        "Process grid ", npy, " x ", npz, " does not match the number of processes ", nprocs
      call MPI_Abort(MPI_COMM_WORLD, 2, ie)
    end if

    !PFFT expects the process mesh in the C order
    call PoisFFT_InitMPIGrid(MPI_COMM_WORLD, int([npz, npy], c_int), cart_comm, ie)
    if (ie /= 0) then
      if (is_root) write(error_unit, '(a)') "PoisFFT_InitMPIGrid failed."
      call MPI_Abort(MPI_COMM_WORLD, 3, ie)
    end if

    call PoisFFT_LocalGridSize(3, int(gnxyz, c_intptr_t), cart_comm, ln, lo, ln2, lo2)
    nxyz = int(ln)
    offs = int(lo)

    if (any(nxyz < 1)) then
      write(error_unit, '(a,i0,a,3(1x,i0))') "Rank ", rank, " has an empty local grid", nxyz
      call MPI_Abort(MPI_COMM_WORLD, 4, ie)
    end if
  end subroutine


  subroutine grid_finish()
    integer :: ie
    call MPI_Finalize(ie)
  end subroutine


  !global relative maximum-norm error over all processes
  real(dp) function global_rel_err(local_err, local_scale) result(res)
    real(dp), intent(in) :: local_err, local_scale
    real(dp) :: loc(3), glob(3)
    integer :: ie

    !NaN is reduced as a separate flag, MPI_MAX of NaN is not portable
    loc = [local_err, local_scale, 0._dp]
    if (ieee_is_nan(local_err)) loc = [0._dp, local_scale, 1._dp]
    call MPI_Allreduce(loc, glob, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ie)
    if (glob(2) == 0) glob(2) = 1
    res = glob(1) / glob(2)
    if (glob(3) > 0) res = ieee_value(1._dp, ieee_quiet_nan)
  end function

end module mpi_test_grid
