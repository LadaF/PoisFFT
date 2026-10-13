#include "precision.h"
!Distributed version of test_eigen_1d.F90 (compiled only with mpi=yes).
!
!  mpirun -n <np> test_mpi_eigen_1d_dp <np> 1 [case]
!
!The 1D solver uses a 1D Cartesian communicator of all processes, the grid
!is split into np contiguous blocks, the first gn mod np blocks are one
!point longer. The first two arguments keep the command line of the 3D
!tests and are not used (the runner passes np and 1). Every process builds
!its part of the global eigenvector and the error is reduced over all
!processes. `case` selects one pair of boundary conditions (1..9), 0 or
!no argument runs all.
program test_mpi_eigen_1d
  use mpi
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: error_unit
  use ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use pfft, only: pfft_create_procmesh_1d
  use PoisFFT, Solver1D => SOLVER1D
  use reference_modes
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)
  real(dp), parameter :: Ls(1) = [2.3_dp]
  integer, parameter :: gn(1) = [7]

  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ds = PoisFFT_DirichletStag, &
                        Ns = PoisFFT_NeumannStag
  integer, parameter :: pairs(2,9) = reshape([P,P, D,D, Ds,Ds, N,N, Ns,Ns, &
                                              N,D, D,N, Ns,Ds, Ds,Ns], [2,9])
  integer, parameter :: approxs(3) = [PoisFFT_Spectral, &
                                      PoisFFT_FiniteDifference2, &
                                      PoisFFT_FiniteDifference4]
  integer :: comm, rank, nprocs, ie, provided, command_case
  integer :: nx(1), offs(1)
  integer :: ip, ia
  character(32) :: arg, grid_name

  call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ie)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ie)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ie)

  command_case = 0
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) command_case
  end if

  ie = pfft_create_procmesh_1d(MPI_COMM_WORLD, int(nprocs, c_int), comm)
  if (ie /= 0) then
    if (rank == 0) write(error_unit, '(a)') "pfft_create_procmesh_1d failed."
    call MPI_Abort(MPI_COMM_WORLD, 3, ie)
  end if

  nx = gn / nprocs
  offs = rank * nx + min(rank, mod(gn, nprocs))
  if (rank < mod(gn(1), nprocs)) nx = nx + 1
  if (any(nx < 1)) then
    write(error_unit, '(a,i0,a)') "Rank ", rank, " has an empty local grid"
    call MPI_Abort(MPI_COMM_WORLD, 4, ie)
  end if

  write(grid_name, '(" np ",i0)') nprocs
  call ut_init("test_mpi_eigen_1d " // PREC_NAME // trim(grid_name), silent=rank/=0)

  do ip = 1, size(pairs, 2)
    if (command_case /= 0 .and. command_case /= ip) cycle
    do ia = 1, size(approxs)
      call ut_section("MPI 1D " // bc_pair_name(pairs(:,ip)) // " " // approx_name(approxs(ia)))
      call test_all_modes(pairs(:,ip), approxs(ia))
    end do
  end do

  call MPI_Finalize(ie)
  call ut_finish()

contains

  subroutine test_all_modes(bcs, approx)
    integer, intent(in) :: bcs(2), approx
    type(Solver1D) :: solver
    real(rp), allocatable :: rhs(:), phi(:), exact(:)
    real(dp) :: vx(gn(1)), lam, err, lam_min
    logical :: ok
    integer :: mx
    character(80) :: label

    allocate(rhs(nx(1)), phi(nx(1)), exact(nx(1)))

    solver = Solver1D(nx, real(Ls, rp), bcs, approx, &
                      gnxyz=gn, offs=offs, mpi_comm=comm, &
                      planning=PoisFFT_Estimate)
    lam_min = min_eigenvalue(bcs, gn, Ls, approx)

    do mx = 1, gn(1)
      call mode_1d(bcs, gn(1), Ls(1), mx, approx, vx, lam, ok)

      rhs = real(vx(offs(1)+1:offs(1)+nx(1)), rp)
      if (lam == 0) then
        exact = 0
      else
        exact = real(vx(offs(1)+1:offs(1)+nx(1)) / lam, rp)
      end if

      call Execute(solver, phi, rhs)

      !normalized as described in reference_modes.f90
      err = global_rel_err(real(maxval(abs(phi - exact)), dp) * lam_min, &
                           real(maxval(abs(rhs)), dp))
      write(label, '("n=",i0," mode ",i0)') gn, mx
      call check_err(trim(label), err, real(tol, dp))
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

end program test_mpi_eigen_1d
