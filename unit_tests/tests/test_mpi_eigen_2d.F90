#include "precision.h"
!Distributed version of test_eigen_2d.F90 (compiled only with mpi=yes).
!
!  mpirun -n <np> test_mpi_eigen_2d_dp <np> 1 [case]
!
!The 2D solver uses a 1D Cartesian communicator of all processes: x is
!never distributed, y is split into np parts. The first two arguments
!keep the command line of the 3D tests and are not used (the runner
!passes np and 1, see the nps option of the SConstruct). The communicator and
!the local block come from PoisFFT_InitMPIGrid and PoisFFT_LocalGridSize for
!the 2D grid, every process builds its part of the
!global eigenvector and the error is reduced over all processes. `case`
!selects one combination of boundary conditions (1..n_combos_2d), 0 or no
!argument runs all.
program test_mpi_eigen_2d
  use mpi
  use iso_c_binding, only: c_int, c_intptr_t
  use iso_fortran_env, only: error_unit
  use ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use PoisFFT, Solver2D => SOLVER2D
  use reference_modes
  use unit_test_utils
  use bc_combos
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)
  real(dp), parameter :: Ls(2) = [2.3_dp, 1.7_dp]
  integer, parameter :: gn(2) = [6, 7]

  integer, parameter :: approxs(3) = [PoisFFT_Spectral, &
                                      PoisFFT_FiniteDifference2, &
                                      PoisFFT_FiniteDifference4]
  integer :: comm, rank, nprocs, ie, provided, command_case
  integer :: nxy(2), offs(2)
  integer(c_intptr_t) :: ln(2), lo(2), ln2(2), lo2(2)
  integer :: ic, ia
  character(32) :: arg, grid_name

  call MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ie)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ie)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ie)

  !the solvers call MPI from inside OpenMP single regions
  !$ if (provided < MPI_THREAD_SERIALIZED) then
  !$   if (rank == 0) write(error_unit, '(a)') "MPI does not provide MPI_THREAD_SERIALIZED."
  !$   call MPI_Abort(MPI_COMM_WORLD, 5, ie)
  !$ end if

  command_case = 0
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) command_case
  end if

  call PoisFFT_InitMPIGrid(MPI_COMM_WORLD, [int(nprocs, c_int)], comm, ie)
  if (ie /= 0) then
    if (rank == 0) write(error_unit, '(a)') "PoisFFT_InitMPIGrid failed."
    call MPI_Abort(MPI_COMM_WORLD, 3, ie)
  end if

  call PoisFFT_LocalGridSize(2, int(gn, c_intptr_t), comm, ln, lo, ln2, lo2)
  nxy = int(ln)
  offs = int(lo)
  if (any(nxy < 1)) then
    write(error_unit, '(a,i0,a,2(1x,i0))') "Rank ", rank, " has an empty local grid", nxy
    call MPI_Abort(MPI_COMM_WORLD, 4, ie)
  end if

  write(grid_name, '(" np ",i0)') nprocs
  call ut_init("test_mpi_eigen_2d " // PREC_NAME // trim(grid_name), silent=rank/=0)

  do ic = 1, n_combos_2d
    if (command_case /= 0 .and. command_case /= ic) cycle
    do ia = 1, size(approxs)
      call ut_section("MPI 2D " // bc_pair_name(combos_2d(:,ic)) // " " // approx_name(approxs(ia)))
      call test_all_modes(combos_2d(:,ic), approxs(ia))
    end do
  end do

  call MPI_Finalize(ie)
  call ut_finish()

contains

  subroutine test_all_modes(bcs, approx)
    integer, intent(in) :: bcs(4), approx
    type(Solver2D) :: solver
    real(rp), allocatable :: rhs(:,:), phi(:,:), exact(:,:)
    real(dp) :: vx(gn(1)), vy(gn(2)), lx, ly, lam
    real(dp) :: err, lam_min
    logical :: ok(2)
    integer :: mx, my, i, j
    character(80) :: label

    allocate(rhs(nxy(1),nxy(2)), phi(nxy(1),nxy(2)), exact(nxy(1),nxy(2)))

    solver = Solver2D(nxy, real(Ls, rp), bcs, approx, &
                      gnxyz=gn, offs=offs, mpi_comm=comm, &
                      planning=PoisFFT_Estimate)
    lam_min = min_eigenvalue(bcs, gn, Ls, approx)

    do my = 1, gn(2)
      call mode_1d(bcs(3:4), gn(2), Ls(2), my, approx, vy, ly, ok(2))
      do mx = 1, gn(1)
        call mode_1d(bcs(1:2), gn(1), Ls(1), mx, approx, vx, lx, ok(1))
        lam = lx + ly

        do j = 1, nxy(2)
          do i = 1, nxy(1)
            rhs(i,j) = real(vx(i+offs(1)) * vy(j+offs(2)), rp)
            if (lam == 0) then
              exact(i,j) = 0
            else
              exact(i,j) = real(vx(i+offs(1)) * vy(j+offs(2)) / lam, rp)
            end if
          end do
        end do

        call Execute(solver, phi, rhs)

        !normalized as described in reference_modes.f90
        err = global_rel_err(real(maxval(abs(phi - exact)), dp) * lam_min, &
                             real(maxval(abs(rhs)), dp))
        write(label, '("n=",i0,"x",i0," mode ",i0,",",i0)') gn, mx, my
        call check_err(trim(label), err, real(tol, dp))
      end do
    end do

    call Finalize(solver)
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

end program test_mpi_eigen_2d
