#include "precision.h"
!Unit test of the error reporting through the optional arguments ierr and
!errmsg of the constructors. PoisFFT_Solver3D_nonuniform_z:
!  - valid arguments give PoisFFT_OK, an empty errmsg and the same solution
!    as a solver constructed without ierr
!  - unsupported boundary conditions give PoisFFT_ErrorBCs
!  - too short z or z_u arrays and non-increasing coordinates give
!    PoisFFT_ErrorGrid
!  - a solver returned with an error can be finalized and a valid solver
!    can be constructed afterwards
!PoisFFT_Solver1D, 2D and 3D: valid arguments as above, unsupported
!boundary conditions, a local grid outside of the global one, a local grid
!different from the global one without MPI and a regular Neumann pair on
!a single point.
!The behaviour without ierr (stop with the message) is tested in
!test_unsupported_bcs.F90, the MPI errors in test_mpi_error_reporting.F90.
program test_error_reporting
  use PoisFFT, Solver_nuz => SOLVER3D_NUZ, Solver1D => SOLVER1D, &
               Solver2D => SOLVER2D, Solver3D => SOLVER3D
  use reference_modes, only: dp, bc_pair_name
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND
  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ns = PoisFFT_NeumannStag, &
                        Ds = PoisFFT_DirichletStag
  integer, parameter :: nxyz(3) = [6, 5, 4]
  real(rp), parameter :: Lxy(2) = [1.3_rp, 0.9_rp]
  real(rp) :: z(nxyz(3)), z_u(0:nxyz(3))
  integer :: k

  call ut_init("test_error_reporting " // PREC_NAME)

  z_u = [(real(k, rp)**1.5_rp / nxyz(3), k = 0, nxyz(3))]
  z = (z_u(1:) + z_u(0:nxyz(3)-1)) / 2

  call ut_section("valid arguments")
  call test_valid([P, P, P, P, Ns, Ns])
  call test_valid([P, P, N, N, Ns, Ns])
  call test_valid([N, N, N, N, Ns, Ns])
  call test_valid([D, D, P, P, Ds, Ds])

  call ut_section("unsupported boundary conditions")
  call test_error([(D, k = 1, 6)], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([P, P, P, P, D, D], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([P, P, P, P, P, D], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  !only staggered conditions in z
  call test_error([P, P, P, P, N, N], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([(N, k = 1, 6)], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([(9, k = 1, 6)], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")

  call ut_section("inconsistent grid")
  call test_error([P, P, P, P, Ns, Ns], z(1:nxyz(3)-1), z_u, PoisFFT_ErrorGrid, "whole grid")
  call test_error([P, P, P, P, Ns, Ns], z, z_u(0:nxyz(3)-1), PoisFFT_ErrorGrid, "whole grid")
  call test_error([P, P, P, P, Ns, Ns], z([1, 3, 2, 4]), z_u, PoisFFT_ErrorGrid, "increasing")
  call test_error([P, P, P, P, Ns, Ns], z, z_u([0, 1, 1, 3, 4]), PoisFFT_ErrorGrid, "increasing")
  !Dirichlet at the face z_u(0) with the first point at the same position
  call test_error([P, P, P, P, Ds, Ns], [z_u(0), z(2:)], z_u, PoisFFT_ErrorGrid, "increasing")

  call ut_section("ierr without errmsg")
  call test_ierr_only()

  call ut_section("PoisFFT_Solver1D/2D/3D valid arguments")
  call test_valid_xd()

  call ut_section("PoisFFT_Solver1D/2D/3D invalid arguments")
  call test_error_1d([P, D], [7], [7], [0], PoisFFT_ErrorBCs, "boundary conditions")
  call test_error_1d([D, D], [0], [0], [0], PoisFFT_ErrorGrid, "inside")
  call test_error_1d([D, D], [7], [7], [1], PoisFFT_ErrorGrid, "inside")
  call test_error_1d([D, D], [7], [8], [0], PoisFFT_ErrorGrid, "whole grid")
  call test_error_1d([N, N], [1], [1], [0], PoisFFT_ErrorGrid, "two points")
  call test_error_2d([P, P, D, Ds], [4, 5], [4, 5], [0, 0], PoisFFT_ErrorBCs, "boundary conditions")
  call test_error_2d([P, P, D, D], [4, 0], [4, 0], [0, 0], PoisFFT_ErrorGrid, "inside")
  call test_error_2d([P, P, D, D], [4, 5], [4, 6], [0, 0], PoisFFT_ErrorGrid, "whole grid")
  call test_error_2d([N, N, P, P], [1, 5], [1, 5], [0, 0], PoisFFT_ErrorGrid, "two points")
  call test_error_3d([D, Ds, P, P, P, P], [4, 5, 3], [4, 5, 3], [0, 0, 0], PoisFFT_ErrorBCs, &
                     "boundary conditions")
  call test_error_3d([P, P, P, P, P, P], [4, 5, 3], [4, 5, 3], [-1, 0, 0], PoisFFT_ErrorGrid, "inside")
  call test_error_3d([P, P, P, P, P, P], [4, 5, 3], [4, 5, 6], [0, 0, 0], PoisFFT_ErrorGrid, "whole grid")
  call test_error_3d([P, P, N, N, D, D], [4, 1, 3], [4, 1, 3], [0, 0, 0], PoisFFT_ErrorGrid, "two points")

  call ut_finish()

contains

  subroutine random_rhs(rhs)
    real(rp), intent(out) :: rhs(:,:,:)
    real(dp), allocatable :: tmp(:,:,:)

    allocate(tmp(size(rhs,1), size(rhs,2), size(rhs,3)))
    call random_number(tmp)
    rhs = real(tmp, rp)
    rhs = rhs - sum(rhs) / size(rhs)
  end subroutine


  subroutine test_valid(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver_nuz) :: s1, s2
    real(rp), dimension(nxyz(1),nxyz(2),nxyz(3)) :: rhs, phi1, phi2
    integer :: ierr
    character(100) :: errmsg

    errmsg = "garbage"
    s1 = Solver_nuz(nxyz, Lxy, z, z_u, bcs, PoisFFT_FiniteDifference2, &
                    planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check(ierr == PoisFFT_OK, bc_pair_name(bcs) // " ierr == PoisFFT_OK")
    call check(len_trim(errmsg) == 0, bc_pair_name(bcs) // " errmsg empty")

    s2 = Solver_nuz(nxyz, Lxy, z, z_u, bcs, PoisFFT_FiniteDifference2, &
                    planning=PoisFFT_Estimate)

    call random_rhs(rhs)
    call Execute(s1, phi1, rhs)
    call Execute(s2, phi2, rhs)
    call check(all(phi1 == phi2), bc_pair_name(bcs) // " same solution as without ierr")

    call Finalize(s1)
    call Finalize(s2)
  end subroutine


  subroutine test_error(bcs, zz, zz_u, expected, msg_part)
    integer, intent(in) :: bcs(6), expected
    real(rp), intent(in) :: zz(:), zz_u(0:)
    character(*), intent(in) :: msg_part
    type(Solver_nuz) :: s
    integer :: ierr
    character(100) :: errmsg
    character(80) :: label

    write(label, '(a," size(z)=",i0," size(z_u)=",i0)') bc_pair_name(bcs), size(zz), size(zz_u)

    s = Solver_nuz(nxyz, Lxy, zz, zz_u, bcs, PoisFFT_FiniteDifference2, &
                   planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check(ierr == expected, trim(label) // " ierr")
    call check(index(errmsg, msg_part) > 0, trim(label) // " errmsg", trim(errmsg))

    !the solver returned with an error can be finalized
    call Finalize(s)

    !and a valid solver can be constructed afterwards
    s = Solver_nuz(nxyz, Lxy, z, z_u, [P, P, P, P, Ns, Ns], PoisFFT_FiniteDifference2, &
                   planning=PoisFFT_Estimate, ierr=ierr)
    call check(ierr == PoisFFT_OK, trim(label) // " valid solver afterwards")
    call Finalize(s)
  end subroutine


  subroutine test_ierr_only()
    type(Solver_nuz) :: s
    integer :: ierr

    s = Solver_nuz(nxyz, Lxy, z, z_u, [(D, k = 1, 6)], ierr=ierr)
    call check(ierr == PoisFFT_ErrorBCs, "unsupported boundary conditions")
    call Finalize(s)
  end subroutine



  !ierr == PoisFFT_OK, errmsg empty and the same solution as without ierr
  subroutine test_valid_xd()
    type(Solver1D) :: a1, b1
    type(Solver2D) :: a2, b2
    type(Solver3D) :: a3, b3
    real(rp) :: r1(7), p1(7), q1(7), r2(4,5), p2(4,5), q2(4,5)
    real(rp) :: r3(4,5,3), p3(4,5,3), q3(4,5,3)
    real(dp) :: t1(7), t2(4,5), t3(4,5,3)
    integer :: ierr
    character(100) :: errmsg

    call random_number(t1); r1 = real(t1, rp)
    call random_number(t2); r2 = real(t2, rp)
    call random_number(t3); r3 = real(t3, rp)

    a1 = Solver1D([7], [1.3_rp], [Ns, Ds], planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check(ierr == PoisFFT_OK .and. len_trim(errmsg) == 0, "1D Ns,Ds ierr and errmsg")
    b1 = Solver1D([7], [1.3_rp], [Ns, Ds], planning=PoisFFT_Estimate)
    call Execute(a1, p1, r1)
    call Execute(b1, q1, r1)
    call check(all(p1 == q1), "1D Ns,Ds same solution as without ierr")
    call Finalize(a1)
    call Finalize(b1)

    a2 = Solver2D([4, 5], Lxy, [P, P, D, N], planning=PoisFFT_Estimate, ierr=ierr, errmsg=errmsg)
    call check(ierr == PoisFFT_OK .and. len_trim(errmsg) == 0, "2D P,P D,N ierr and errmsg")
    b2 = Solver2D([4, 5], Lxy, [P, P, D, N], planning=PoisFFT_Estimate)
    call Execute(a2, p2, r2)
    call Execute(b2, q2, r2)
    call check(all(p2 == q2), "2D P,P D,N same solution as without ierr")
    call Finalize(a2)
    call Finalize(b2)

    a3 = Solver3D([4, 5, 3], [Lxy, 0.7_rp], [P, P, N, N, Ds, Ds], &
                  gnxyz=[4, 5, 3], offs=[0, 0, 0], planning=PoisFFT_Estimate, &
                  ierr=ierr, errmsg=errmsg)
    call check(ierr == PoisFFT_OK .and. len_trim(errmsg) == 0, "3D P,P N,N Ds,Ds ierr and errmsg")
    b3 = Solver3D([4, 5, 3], [Lxy, 0.7_rp], [P, P, N, N, Ds, Ds], planning=PoisFFT_Estimate)
    call Execute(a3, p3, r3)
    call Execute(b3, q3, r3)
    call check(all(p3 == q3), "3D P,P N,N Ds,Ds same solution as without ierr")
    call Finalize(a3)
    call Finalize(b3)
  end subroutine


  subroutine check_xd_error(label, ierr, errmsg, expected, msg_part)
    character(*), intent(in) :: label, errmsg, msg_part
    integer, intent(in) :: ierr, expected

    call check(ierr == expected, label // " ierr")
    call check(index(errmsg, msg_part) > 0, label // " errmsg", trim(errmsg))
  end subroutine


  subroutine test_error_1d(bcs, n, gn, off, expected, msg_part)
    integer, intent(in) :: bcs(2), n(1), gn(1), off(1), expected
    character(*), intent(in) :: msg_part
    type(Solver1D) :: s
    integer :: ierr
    character(100) :: errmsg
    character(80) :: label

    write(label, '("1D ",a," n=",i0," gn=",i0," off=",i0)') bc_pair_name(bcs), n, gn, off
    s = Solver1D(n, [1.3_rp], bcs, gnxyz=gn, offs=off, planning=PoisFFT_Estimate, &
                 ierr=ierr, errmsg=errmsg)
    call check_xd_error(trim(label), ierr, errmsg, expected, msg_part)
    !the solver returned with an error can be finalized
    call Finalize(s)
  end subroutine


  subroutine test_error_2d(bcs, n, gn, off, expected, msg_part)
    integer, intent(in) :: bcs(4), n(2), gn(2), off(2), expected
    character(*), intent(in) :: msg_part
    type(Solver2D) :: s
    integer :: ierr
    character(100) :: errmsg
    character(80) :: label

    write(label, '("2D ",a," n=",i0,"x",i0," gn=",i0,"x",i0)') bc_pair_name(bcs), n, gn
    s = Solver2D(n, Lxy, bcs, gnxyz=gn, offs=off, planning=PoisFFT_Estimate, &
                 ierr=ierr, errmsg=errmsg)
    call check_xd_error(trim(label), ierr, errmsg, expected, msg_part)
    call Finalize(s)
  end subroutine


  subroutine test_error_3d(bcs, n, gn, off, expected, msg_part)
    integer, intent(in) :: bcs(6), n(3), gn(3), off(3), expected
    character(*), intent(in) :: msg_part
    type(Solver3D) :: s
    integer :: ierr
    character(100) :: errmsg
    character(80) :: label

    write(label, '("3D ",a," n=",i0,2("x",i0)," gn=",i0,2("x",i0)," off=",i0,2(",",i0))') &
          bc_pair_name(bcs), n, gn, off
    s = Solver3D(n, [Lxy, 0.7_rp], bcs, gnxyz=gn, offs=off, planning=PoisFFT_Estimate, &
                 ierr=ierr, errmsg=errmsg)
    call check_xd_error(trim(label), ierr, errmsg, expected, msg_part)
    call Finalize(s)
  end subroutine

end program test_error_reporting
