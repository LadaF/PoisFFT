!Unit test of the rejection of unsupported boundary conditions and of
!inconsistent grid arguments when no ierr argument is given.
!
!  test_unsupported_bcs <case>
!
!Case 1..n_cases constructs one solver with an unsupported combination of
!boundary conditions (cases 18, 19: inconsistent z grid of the non-uniform
!solver, 20, 21: inconsistent grid of the 3D and 1D solvers). The library must stop with the error message, e.g.
!"Unknown combination of boundary conditions."; if the constructor returns,
!the program prints NOT REJECTED and fails. Case 0 constructs valid solvers
!of every kind and must finish normally (prints "all accepted").
!The test runner checks the output of every case, `stop` with a message
!does not give a non-zero exit code.
program test_unsupported_bcs
  use PoisFFT, Solver1D => PoisFFT_Solver1D_DP, &
               Solver2D => PoisFFT_Solver2D_DP, &
               Solver3D => PoisFFT_Solver3D_DP, &
               Solver_nuz => PoisFFT_Solver3D_nonuniform_z_DP
  use iso_fortran_env, only: dp => real64, output_unit
  implicit none

  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ds = PoisFFT_DirichletStag, &
                        Ns = PoisFFT_NeumannStag
  integer, parameter :: bad = 9   !not a boundary condition at all

  integer, parameter :: n_cases = 21
  integer :: icase, k
  character(16) :: arg
  real(dp) :: z(4), z_u(0:4)

  z_u = [(0.25_dp * k, k = 0, 4)]
  z = (z_u(1:4) + z_u(0:3)) / 2

  icase = 0
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) icase
  end if

  select case (icase)
  case (0)
    call try_1d([P, P]);  call try_1d([D, N]);  call try_1d([Ns, Ds])
    call try_2d([(Ds, k = 1, 4)]);      call try_2d([P, P, D, N])
    call try_2d([D, D, Ns, Ns])
    call try_3d([P, P, P, P, Ds, Ns]);  call try_3d([D, D, D, D, N, N])
    call try_3d([D, D, P, P, P, P]);    call try_3d([N, D, Ns, Ns, P, P])
    call try_nuz([P, P, P, P, Ns, Ns]); call try_nuz([D, D, P, P, Ds, Ns])
    write(output_unit, '(a)') "all accepted"
    stop
  !1D: periodic with non-periodic, mixed regular and staggered, invalid code
  case (1);  call try_1d([P, D])
  case (2);  call try_1d([Ns, P])
  case (3);  call try_1d([D, Ds])
  case (4);  call try_1d([N, Ns])
  case (5);  call try_1d([bad, bad])
  !2D: an invalid pair in one direction
  case (6);  call try_2d([P, P, D, Ds])
  case (7);  call try_2d([D, P, N, N])
  case (8);  call try_2d([Ds, Ds, Ns, N])
  case (9);  call try_2d([P, D, P, P])
  !3D: an invalid pair in one direction
  case (10); call try_3d([D, D, N, Ns, P, P])
  case (11); call try_3d([P, P, P, P, D, Ds])
  case (12); call try_3d([P, D, P, P, P, P])
  case (13); call try_3d([N, N, D, D, P, N])
  case (14); call try_3d([(bad, k = 1, 6)])
  !non-uniform z: only staggered conditions in z, any supported pairs in x and y
  case (15); call try_nuz([P, P, P, P, N, N])
  case (16); call try_nuz([P, P, P, P, D, D])
  case (17); call try_nuz([P, P, P, P, P, D])
  !non-uniform z: inconsistent grid
  case (18); call try_nuz_grid(z(1:3), z_u)
  case (19); call try_nuz_grid(z([1, 3, 2, 4]), z_u)
  !3D: the regular Neumann pair on a single point; 1D: a part of the grid
  !without MPI
  case (20); call try_3d_n([P, P, N, N, D, D], [8, 1, 4])
  case (21); call try_1d_part()
  case default
    write(output_unit, '(a,i0)') "unknown test case ", icase
    error stop 2
  end select

  !only reached when the solver accepted the boundary conditions
  write(output_unit, '(a,i0)') "NOT REJECTED, case ", icase
  error stop 1

contains

  subroutine try_1d(bcs)
    integer, intent(in) :: bcs(2)
    type(Solver1D) :: s
    s = Solver1D([8], [1._dp], bcs, planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

  subroutine try_2d(bcs)
    integer, intent(in) :: bcs(4)
    type(Solver2D) :: s
    s = Solver2D([8, 6], [1._dp, 1._dp], bcs, planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

  subroutine try_3d(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver3D) :: s
    s = Solver3D([8, 6, 4], [1._dp, 1._dp, 1._dp], bcs, planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

  subroutine try_3d_n(bcs, n)
    integer, intent(in) :: bcs(6), n(3)
    type(Solver3D) :: s
    s = Solver3D(n, [1._dp, 1._dp, 1._dp], bcs, planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

  subroutine try_1d_part()
    type(Solver1D) :: s
    s = Solver1D([4], [1._dp], [D, D], gnxyz=[8], offs=[2], planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

  subroutine try_nuz(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver_nuz) :: s
    s = Solver_nuz([8, 6, 4], [1._dp, 1._dp], z, z_u, bcs, planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

  subroutine try_nuz_grid(zz, zz_u)
    real(dp), intent(in) :: zz(:), zz_u(0:)
    type(Solver_nuz) :: s
    s = Solver_nuz([8, 6, 4], [1._dp, 1._dp], zz, zz_u, [P, P, P, P, Ns, Ns], &
                   planning=PoisFFT_Estimate)
    call Finalize(s)
  end subroutine

end program test_unsupported_bcs
