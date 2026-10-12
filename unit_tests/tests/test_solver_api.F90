#include "precision.h"
!Unit test of the general behaviour of the solver interface
!(PoisFFT_SolverXD constructors, Execute, Finalize) independent of the
!particular eigenvalues:
!  - ghost layers of any width around Phi and RHS
!  - RHS is not modified by Execute
!  - repeated calls of Execute with the same solver
!  - bitwise reproducibility of two solvers created with PoisFFT_Estimate
!  - PoisFFT_Measure plans give the same result up to rounding
!  - linearity of the solution operator
!  - zero mean of the solution of the singular periodic and Neumann problems
!  - Finalize and creation of a new solver of a different size
program test_solver_api
  use PoisFFT, Solver1D => SOLVER1D, &
               Solver2D => SOLVER2D, &
               Solver3D => SOLVER3D
  use reference_modes, only: dp, bc_pair_name
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)

  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ds = PoisFFT_DirichletStag, &
                        Ns = PoisFFT_NeumannStag
  integer, parameter :: ns3(3) = [6, 5, 7]
  real(rp), parameter :: Ls3(3) = [1.3_rp, 2.1_rp, 0.7_rp]
  integer, parameter :: n_combos = 8
  integer, parameter :: combos(6, n_combos) = reshape([ &
      P, P, P, P, P, P,  Ns,Ns,Ns,Ns,Ns,Ns,  D, D, D, D, D, D, &
      P, P, P, P, Ns,Ns, P, P, Ns,Ns,Ns,Ns,  Ds,Ds,Ds,Ds,P, P, &
      Ns,Ns,P, P, Ns,Ns, D, D, P, P, N, N], [6, n_combos])
  integer :: ic

  call ut_init("test_solver_api " // PREC_NAME)

  do ic = 1, n_combos
    call ut_section("3D " // bc_pair_name(combos(:,ic)))
    call test_ghost_layers_3d(combos(:,ic))
    call test_repeated_and_reproducible_3d(combos(:,ic))
    call test_linearity_3d(combos(:,ic))
  end do

  call ut_section("3D singular problems")
  call test_zero_mean_3d([(P, ic = 1, 6)])
  call test_zero_mean_3d([(Ns, ic = 1, 6)])
  call test_zero_mean_3d([P, P, P, P, Ns, Ns])
  call test_zero_mean_3d([P, P, Ns, Ns, Ns, Ns])
  call test_zero_mean_3d([Ns, Ns, P, P, Ns, Ns])

  call ut_section("2D")
  call test_ghost_layers_2d([(P, ic = 1, 4)])
  call test_ghost_layers_2d([(Ds, ic = 1, 4)])
  call test_ghost_layers_2d([(N, ic = 1, 4)])

  call ut_section("1D")
  call test_ghost_layers_1d([P, P])
  call test_ghost_layers_1d([D, N])
  call test_ghost_layers_1d([Ns, Ns])

  call ut_section("Finalize and re-create")
  call test_recreate()

  call ut_finish()

contains

  subroutine random_rhs_3d(rhs, bcs)
    real(rp), intent(out) :: rhs(:,:,:)
    integer, intent(in) :: bcs(6)
    real(dp), allocatable :: tmp(:,:,:)

    allocate(tmp(size(rhs,1), size(rhs,2), size(rhs,3)))
    call random_number(tmp)
    rhs = real(tmp - 0.5_dp, rp)
    !compatibility condition of the singular problems
    if (all(bcs == P .or. bcs == Ns)) rhs = rhs - sum(rhs) / size(rhs)
  end subroutine


  subroutine test_ghost_layers_3d(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver3D) :: solver
    real(rp) :: rhs(ns3(1),ns3(2),ns3(3)), phi(ns3(1),ns3(2),ns3(3))
    real(rp), allocatable :: rhs_g(:,:,:), phi_g(:,:,:)
    integer :: gr, gp
    character(40) :: label

    solver = Solver3D(ns3, Ls3, bcs, PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)

    call random_rhs_3d(rhs, bcs)
    call Execute(solver, phi, rhs)

    do gr = 0, 2
      do gp = 0, 3
        if (gr == 0 .and. gp == 0) cycle
        allocate(rhs_g(1-gr:ns3(1)+gr, 1-gr:ns3(2)+gr, 1-gr:ns3(3)+gr))
        allocate(phi_g(1-gp:ns3(1)+gp, 1-gp:ns3(2)+gp, 1-gp:ns3(3)+gp))
        rhs_g = 1e3_rp   !garbage in the ghost points must be ignored
        rhs_g(1:ns3(1), 1:ns3(2), 1:ns3(3)) = rhs
        phi_g = 0

        call Execute(solver, phi_g, rhs_g)

        write(label, '("ghosts RHS=",i0," Phi=",i0)') gr, gp
        call check_err(trim(label), rel_err(phi_g(1:ns3(1), 1:ns3(2), 1:ns3(3)), phi), tol)
        deallocate(rhs_g, phi_g)
      end do
    end do

    call Finalize(solver)
  end subroutine


  subroutine test_repeated_and_reproducible_3d(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver3D) :: s1, s2, s3
    real(rp), dimension(ns3(1),ns3(2),ns3(3)) :: rhs, rhs_copy, phi1, phi2, phi3, phi_m

    s1 = Solver3D(ns3, Ls3, bcs, PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)
    s2 = Solver3D(ns3, Ls3, bcs, PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)
    s3 = Solver3D(ns3, Ls3, bcs, PoisFFT_FiniteDifference2)

    call random_rhs_3d(rhs, bcs)
    rhs_copy = rhs

    call Execute(s1, phi1, rhs)
    call check(all(rhs == rhs_copy), "RHS not modified by Execute")

    call Execute(s1, phi2, rhs)
    call check(all(phi1 == phi2), "repeated Execute gives identical result")

    call Execute(s2, phi2, rhs)
    call check(all(phi1 == phi2), "PoisFFT_Estimate solvers are bitwise reproducible")

    call Execute(s3, phi3, rhs)
    call check_err("PoisFFT_Measure solver agrees", rel_err(phi3, phi1), tol)

    !a different RHS after the first one, the solver keeps no state
    call random_rhs_3d(rhs, bcs)
    call Execute(s1, phi_m, rhs)
    call Execute(s2, phi2, rhs)
    call check(all(phi_m == phi2), "second RHS, state not kept between calls")

    call Finalize(s1)
    call Finalize(s2)
    call Finalize(s3)
  end subroutine


  subroutine test_linearity_3d(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver3D) :: solver
    real(rp), dimension(ns3(1),ns3(2),ns3(3)) :: f, g, phi_f, phi_g, phi_fg
    real(rp), parameter :: a = 2.5_rp, b = -0.75_rp

    solver = Solver3D(ns3, Ls3, bcs, PoisFFT_Spectral, planning=PoisFFT_Estimate)

    call random_rhs_3d(f, bcs)
    call random_rhs_3d(g, bcs)
    call Execute(solver, phi_f, f)
    call Execute(solver, phi_g, g)
    call Execute(solver, phi_fg, a*f + b*g)

    call check_err("linearity", rel_err(phi_fg, a*phi_f + b*phi_g), tol)

    call Finalize(solver)
  end subroutine


  subroutine test_zero_mean_3d(bcs)
    integer, intent(in) :: bcs(6)
    type(Solver3D) :: solver
    real(rp), dimension(ns3(1),ns3(2),ns3(3)) :: rhs, phi

    solver = Solver3D(ns3, Ls3, bcs, PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)

    call random_rhs_3d(rhs, bcs)
    !a constant added to the RHS must be ignored (projected out)
    rhs = rhs + 3
    call Execute(solver, phi, rhs)

    call check_err(bc_pair_name(bcs) // " solution has zero mean", &
                   abs(sum(phi)) / (size(phi) * maxval(abs(phi))), tol)

    call Finalize(solver)
  end subroutine


  subroutine test_ghost_layers_2d(bcs)
    integer, intent(in) :: bcs(4)
    integer, parameter :: nxyz(2) = [7, 6]
    type(Solver2D) :: solver
    real(rp) :: rhs(nxyz(1),nxyz(2)), phi(nxyz(1),nxyz(2))
    real(rp) :: rhs_g(-1:nxyz(1)+2, -1:nxyz(2)+2), phi_g(0:nxyz(1)+1, 0:nxyz(2)+1)
    real(dp) :: tmp(nxyz(1),nxyz(2))

    solver = Solver2D(nxyz, [1.1_rp, 0.9_rp], bcs, PoisFFT_FiniteDifference2, &
                      planning=PoisFFT_Estimate)

    call random_number(tmp)
    rhs = real(tmp, rp) - sum(real(tmp, rp)) / size(tmp)
    call Execute(solver, phi, rhs)

    rhs_g = 1e3_rp
    rhs_g(1:nxyz(1), 1:nxyz(2)) = rhs
    call Execute(solver, phi_g, rhs_g)

    call check_err(bc_pair_name(bcs) // " ghost layers", rel_err(phi_g(1:nxyz(1), 1:nxyz(2)), phi), tol)

    call Finalize(solver)
  end subroutine


  subroutine test_ghost_layers_1d(bcs)
    integer, intent(in) :: bcs(2)
    integer, parameter :: n = 9
    type(Solver1D) :: solver
    real(rp) :: rhs(n), phi(n), rhs_g(-2:n+3), phi_g(0:n+1)
    real(dp) :: tmp(n)

    solver = Solver1D([n], [1.1_rp], bcs, PoisFFT_FiniteDifference2, &
                      planning=PoisFFT_Estimate)

    call random_number(tmp)
    rhs = real(tmp, rp) - sum(real(tmp, rp)) / n
    call Execute(solver, phi, rhs)

    rhs_g = 1e3_rp
    rhs_g(1:n) = rhs
    call Execute(solver, phi_g, rhs_g)

    call check_err(bc_pair_name(bcs) // " ghost layers", rel_err(phi_g(1:n), phi), tol)

    call Finalize(solver)
  end subroutine


  subroutine test_recreate()
    integer, parameter :: bcs(6) = P
    type(Solver3D) :: solver
    real(rp), allocatable :: rhs(:,:,:), phi(:,:,:), phi_ref(:,:,:)
    integer :: i, n

    do i = 1, 3
      n = 3 + 2*i
      allocate(rhs(n, n+1, n+2), phi(n, n+1, n+2), phi_ref(n, n+1, n+2))
      call random_rhs_3d(rhs, bcs)

      solver = Solver3D([n, n+1, n+2], Ls3, bcs, &
                        PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)
      call Execute(solver, phi_ref, rhs)
      call Finalize(solver)

      !the same problem again after Finalize
      solver = Solver3D([n, n+1, n+2], Ls3, bcs, &
                        PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)
      call Execute(solver, phi, rhs)
      call Finalize(solver)

      call check(all(phi == phi_ref), "re-created solver gives identical result")
      deallocate(rhs, phi, phi_ref)
    end do
  end subroutine

end program test_solver_api
