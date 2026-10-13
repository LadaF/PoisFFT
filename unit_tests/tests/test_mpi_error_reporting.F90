#include "precision.h"
!Unit test of the MPI errors reported through ierr and errmsg by the
!constructors (compiled only with mpi=yes).
!
!  mpirun -n <npy*npz> test_mpi_error_reporting_dp <npy> <npz> [case]
!
!  - valid arguments give PoisFFT_OK on all processes
!  - a missing mpi_comm and a communicator that is not Cartesian give
!    PoisFFT_ErrorMPI
!  - a local block that differs from the PFFT decomposition on the first
!    process only gives PoisFFT_ErrorGrid on all processes, with the
!    message about another process on the others, and no process hangs
!  - a 1D grid outside of the global one gives PoisFFT_ErrorGrid
!Every check must hold on all processes. The case argument is not used.
program test_mpi_error_reporting
  use mpi
  use iso_c_binding, only: c_int, c_intptr_t
  use PoisFFT, Solver1D => SOLVER1D, Solver2D => SOLVER2D, &
               Solver3D => SOLVER3D, Solver_nuz => SOLVER3D_NUZ
  use reference_modes, only: dp
  use unit_test_utils
  use mpi_test_grid
  implicit none

  integer, parameter :: rp = RP_KIND
  integer, parameter :: gn(3) = [5, 6, 7]
  real(rp), parameter :: Ls(3) = [2.3_rp, 1.7_rp, 1.1_rp]
  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        Ns = PoisFFT_NeumannStag
  character(32) :: grid_name

  call grid_init(gn)

  write(grid_name, '(" grid ",i0,"x",i0)') npy, npz
  call ut_init("test_mpi_error_reporting " // PREC_NAME // trim(grid_name), silent=.not.is_root)

  call ut_section("3D")
  call test_3d()

  call ut_section("non-uniform z")
  call test_nuz()

  call ut_section("2D and 1D")
  call test_2d_1d()

  call grid_finish()
  call ut_finish()

contains

  !the condition on all processes
  subroutine check_all(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    logical :: glob
    integer :: ie

    call MPI_Allreduce(condition, glob, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ie)
    call check(glob, label)
  end subroutine


  subroutine test_3d()
    type(Solver3D) :: s
    integer :: ierr, n(3)
    character(100) :: errmsg

    s = Solver3D(nxyz, Ls, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, mpi_comm=cart_comm, &
                 planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_OK .and. len_trim(errmsg) == 0, "valid arguments")
    call Finalize(s)

    s = Solver3D(nxyz, Ls, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, &
                 planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_ErrorMPI .and. index(errmsg, "required") > 0, &
                   "missing mpi_comm")
    call Finalize(s)

    s = Solver3D(nxyz, Ls, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, mpi_comm=MPI_COMM_WORLD, &
                 planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_ErrorMPI .and. index(errmsg, "Cartesian") > 0, &
                   "communicator not Cartesian")
    call Finalize(s)

    !x is never distributed, the first process takes one point less
    n = nxyz
    if (rank == 0) n(1) = n(1) - 1
    s = Solver3D(n, Ls, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, mpi_comm=cart_comm, &
                 planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_ErrorGrid, "wrong block on the first process: ierr")
    if (rank == 0) then
      call check_all(index(errmsg, "PFFT decomposition") > 0, &
                     "wrong block on the first process: errmsg")
    else
      call check_all(index(errmsg, "another process") > 0, &
                     "wrong block on the first process: errmsg")
    end if
    call Finalize(s)

    !a valid solver afterwards
    s = Solver3D(nxyz, Ls, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, mpi_comm=cart_comm, &
                 planning=PoisFFT_Estimate, ierr=ierr)
    call check_all(ierr == PoisFFT_OK, "valid solver afterwards")
    call Finalize(s)
  end subroutine


  subroutine test_nuz()
    type(Solver_nuz) :: s
    real(rp) :: z(gn(3)), z_u(0:gn(3))
    integer :: ierr, k
    character(100) :: errmsg

    z_u = [(Ls(3) * k / gn(3), k = 0, gn(3))]
    z = (z_u(1:) + z_u(0:gn(3)-1)) / 2

    s = Solver_nuz(nxyz, Ls(1:2), z, z_u, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, &
                   mpi_comm=cart_comm, planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_OK .and. len_trim(errmsg) == 0, "valid arguments")
    call Finalize(s)

    s = Solver_nuz(nxyz, Ls(1:2), z, z_u, [P, P, D, D, Ns, Ns], gnxyz=gn, offs=offs, &
                   planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_ErrorMPI .and. index(errmsg, "required") > 0, &
                   "missing mpi_comm")
    call Finalize(s)
  end subroutine


  subroutine test_2d_1d()
    type(Solver2D) :: s2
    type(Solver1D) :: s1
    integer :: comm1d, ierr, ie
    integer(c_intptr_t) :: ln(2), lo(2), ln2(2), lo2(2)
    character(100) :: errmsg
    !y = 7 gives non-empty PFFT blocks on up to 4 processes
    integer, parameter :: gn2(2) = [gn(1), gn(3)]

    call PoisFFT_InitMPIGrid(MPI_COMM_WORLD, [int(nprocs, c_int)], comm1d, ie)
    call PoisFFT_LocalGridSize(2, int(gn2, c_intptr_t), comm1d, ln, lo, ln2, lo2)

    s2 = Solver2D(int(ln), Ls(1:2), [P, P, D, D], gnxyz=gn2, offs=int(lo), &
                  mpi_comm=comm1d, planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_OK .and. len_trim(errmsg) == 0, "2D valid arguments")
    call Finalize(s2)

    s2 = Solver2D(int(ln), Ls(1:2), [P, P, D, D], gnxyz=gn2, offs=int(lo), &
                  mpi_comm=MPI_COMM_WORLD, planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_ErrorMPI .and. index(errmsg, "Cartesian") > 0, &
                   "2D communicator not Cartesian")
    call Finalize(s2)

    !one point on every process, the last one beyond the global grid
    s1 = Solver1D([1], Ls(1:1), [D, D], gnxyz=[nprocs - 1], offs=[rank], &
                  mpi_comm=comm1d, planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    if (rank == nprocs - 1) then
      call check_all(ierr == PoisFFT_ErrorGrid .and. index(errmsg, "inside") > 0, &
                     "1D grid outside of the global one")
    else
      call check_all(ierr == PoisFFT_ErrorGrid .and. index(errmsg, "another process") > 0, &
                     "1D grid outside of the global one")
    end if
    call Finalize(s1)

    s1 = Solver1D([1], Ls(1:1), [D, D], gnxyz=[nprocs], offs=[rank], &
                  mpi_comm=MPI_COMM_WORLD, planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check_all(ierr == PoisFFT_ErrorMPI .and. index(errmsg, "Cartesian") > 0, &
                   "1D communicator not Cartesian")
    call Finalize(s1)

    call MPI_Comm_free(comm1d, ie)
  end subroutine

end program test_mpi_error_reporting
