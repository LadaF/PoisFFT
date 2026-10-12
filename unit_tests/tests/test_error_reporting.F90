#include "precision.h"
!Unit test of the error reporting through the optional arguments ierr and
!errmsg of the PoisFFT_Solver3D_nonuniform_z constructor:
!  - valid arguments give PoisFFT_OK, an empty errmsg and the same solution
!    as a solver constructed without ierr
!  - unsupported boundary conditions give PoisFFT_ErrorBCs
!  - too short z or z_u arrays and non-increasing coordinates give
!    PoisFFT_ErrorGrid
!  - a solver returned with an error can be finalized and a valid solver
!    can be constructed afterwards
!The behaviour without ierr (stop with the message) is tested in
!test_unsupported_bcs.F90.
program test_error_reporting
  use PoisFFT, Solver_nuz => SOLVER3D_NUZ
  use reference_modes, only: dp, bc_pair_name
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND
  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ns = PoisFFT_NeumannStag
  integer, parameter :: nxyz(3) = [6, 5, 4]
  real(rp), parameter :: Lxy(2) = [1.3_rp, 0.9_rp]
  real(rp) :: z(nxyz(3)), z_u(0:nxyz(3))
  integer :: k

  call ut_init("test_error_reporting " // PREC_NAME)

  z_u = [(real(k, rp)**1.5_rp / nxyz(3), k = 0, nxyz(3))]
  z = (z_u(1:) + z_u(0:nxyz(3)-1)) / 2

  call ut_section("valid arguments")
  call test_valid([P, P, P, P, Ns, Ns])
  call test_valid([(N, k = 1, 6)])

  call ut_section("unsupported boundary conditions")
  call test_error([(D, k = 1, 6)], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([P, P, P, P, D, D], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([P, P, P, P, P, D], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([P, P, N, N, Ns, Ns], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")
  call test_error([(9, k = 1, 6)], z, z_u, PoisFFT_ErrorBCs, "boundary conditions")

  call ut_section("inconsistent grid")
  call test_error([P, P, P, P, Ns, Ns], z(1:nxyz(3)-1), z_u, PoisFFT_ErrorGrid, "whole grid")
  call test_error([P, P, P, P, Ns, Ns], z, z_u(0:nxyz(3)-1), PoisFFT_ErrorGrid, "whole grid")
  call test_error([P, P, P, P, Ns, Ns], z([1, 3, 2, 4]), z_u, PoisFFT_ErrorGrid, "increasing")
  call test_error([P, P, P, P, Ns, Ns], z, z_u([0, 1, 1, 3, 4]), PoisFFT_ErrorGrid, "increasing")

  call ut_section("ierr without errmsg")
  call test_ierr_only()

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

end program test_error_reporting
