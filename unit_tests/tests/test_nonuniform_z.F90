#include "precision.h"
!Unit test of PoisFFT_Solver3D_nonuniform_z (FFT in x and y, tridiagonal
!solver in z on a non-uniform staggered grid).
!
!The exact solution is phi(i,j,k) = vx(i) vy(j) f(k), where vx, vy are
!1D eigenvectors (reference_modes) and f is an arbitrary profile.
!The right-hand side is vx vy (A f + (lx + ly) f) with A the 2nd order
!finite difference operator on the non-uniform grid with zero flux at the
!boundaries:
!  (A f)_k = [ (f_k+1 - f_k) / (z_k+1 - z_k)
!            - (f_k - f_k-1) / (z_k - z_k-1) ] / (z_u_k - z_u_k-1)
!For the x-y null mode the solution is unique only up to a constant,
!the comparison is then done after subtracting the mean.
!It also checks that on a uniform grid the solution agrees with
!PoisFFT_Solver3D using the FD2 approximation in z.
program test_nonuniform_z
  use PoisFFT, Solver3D => SOLVER3D, Solver_nuz => SOLVER3D_NUZ
  use reference_modes
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 1000 * epsilon(1._rp)

  integer, parameter :: P = PoisFFT_Periodic, N = PoisFFT_Neumann, &
                        Ns = PoisFFT_NeumannStag
  real(dp), parameter :: Lxy(2) = [2.3_dp, 1.7_dp], Lz = 1.1_dp
  !only the z boundary condition type is irrelevant here, the z operator
  !is always given by z and z_u
  integer, parameter :: n_combos = 4
  integer, parameter :: combos(6, n_combos) = reshape([ &
      P, P, P, P, Ns,Ns, &
      P, P, P, P, N, N,  &
      Ns,Ns,Ns,Ns,Ns,Ns, &
      N, N, N, N, N, N], [6, n_combos])
  integer, parameter :: approxs(2) = [PoisFFT_Spectral, PoisFFT_FiniteDifference2]
  integer :: ic, ia

  call ut_init("test_nonuniform_z " // PREC_NAME)

  do ic = 1, n_combos
    do ia = 1, size(approxs)
      call ut_section("non-uniform z " // bc_pair_name(combos(:,ic)) // " " // &
                      approx_name(approxs(ia)))
      call test_modes(combos(:,ic), approxs(ia), [4, 5, 7], stretched=.true.)
      call test_modes(combos(:,ic), approxs(ia), [5, 4, 1], stretched=.true.)
      call test_modes(combos(:,ic), approxs(ia), [3, 4, 6], stretched=.false.)
    end do
  end do

  call ut_section("uniform z grid against PoisFFT_Solver3D")
  call test_uniform_against_3d([P, P, P, P, Ns, Ns])
  call test_uniform_against_3d([(Ns, ic = 1, 6)])

  call ut_finish()

contains

  subroutine make_grid(nz, stretched, z, z_u)
    integer, intent(in) :: nz
    logical, intent(in) :: stretched
    real(dp), intent(out) :: z(1:nz), z_u(0:nz)
    integer :: k

    do k = 0, nz
      if (stretched) then
        z_u(k) = Lz * (real(k, dp) / nz)**1.5_dp
      else
        z_u(k) = Lz * k / nz
      end if
    end do
    z = (z_u(1:nz) + z_u(0:nz-1)) / 2
  end subroutine


  function apply_z(f, z, z_u) result(res)
    real(dp), intent(in) :: f(:), z(:), z_u(0:)
    real(dp) :: res(size(f))
    integer :: k, nz

    nz = size(f)
    res = 0
    do k = 1, nz
      if (k < nz) res(k) = res(k) + (f(k+1) - f(k)) / (z(k+1) - z(k))
      if (k > 1)  res(k) = res(k) - (f(k) - f(k-1)) / (z(k) - z(k-1))
      res(k) = res(k) / (z_u(k) - z_u(k-1))
    end do
  end function


  subroutine test_modes(bcs, approx, nxyz, stretched)
    integer, intent(in) :: bcs(6), approx, nxyz(3)
    logical, intent(in) :: stretched
    type(Solver_nuz) :: solver
    real(dp) :: z(nxyz(3)), z_u(0:nxyz(3)), f(nxyz(3)), af(nxyz(3))
    real(dp) :: vx(nxyz(1)), vy(nxyz(2)), lx, ly
    real(rp), allocatable :: rhs(:,:,:), phi(:,:,:), exact(:,:,:)
    logical :: ok
    integer :: mx, my, j, k
    character(80) :: label

    allocate(rhs(nxyz(1),nxyz(2),nxyz(3)), phi(nxyz(1),nxyz(2),nxyz(3)), exact(nxyz(1),nxyz(2),nxyz(3)))

    call make_grid(nxyz(3), stretched, z, z_u)

    solver = Solver_nuz(nxyz, real(Lxy, rp), real(z, rp), real(z_u, rp), bcs, approx, &
                        planning=PoisFFT_Estimate)

    call random_number(f)
    af = apply_z(f, z, z_u)

    do my = 1, nxyz(2)
      call mode_1d(bcs(3:4), nxyz(2), Lxy(2), my, approx, vy, ly, ok)
      do mx = 1, nxyz(1)
        call mode_1d(bcs(1:2), nxyz(1), Lxy(1), mx, approx, vx, lx, ok)

        do k = 1, nxyz(3)
          do j = 1, nxyz(2)
            rhs(:,j,k) = real(vx * vy(j) * (af(k) + (lx + ly) * f(k)), rp)
            exact(:,j,k) = real(vx * vy(j) * f(k), rp)
          end do
        end do

        call Execute(solver, phi, rhs)

        if (lx + ly == 0) then
          exact = exact - sum(exact) / size(exact)
          phi = phi - sum(phi) / size(phi)
        end if

        write(label, '("n=",i0,"x",i0,"x",i0,a," mode ",i0,",",i0)') nxyz, &
              merge(" stretched", " uniform  ", stretched), mx, my
        !the x-y null mode is a singular tridiagonal system in z made regular
        !by modifying its first row, it is worse conditioned on stretched grids
        call check_err(trim(label), rel_err(phi, exact), &
                       merge(10*tol, tol, lx + ly == 0))
      end do
    end do

    call Finalize(solver)
  end subroutine


  subroutine test_uniform_against_3d(bcs)
    integer, intent(in) :: bcs(6)
    integer, parameter :: nxyz(3) = [6, 5, 7]
    type(Solver_nuz) :: s_nuz
    type(Solver3D) :: s_3d
    real(dp) :: z(nxyz(3)), z_u(0:nxyz(3)), tmp(nxyz(1),nxyz(2),nxyz(3))
    real(rp), dimension(nxyz(1),nxyz(2),nxyz(3)) :: rhs, phi_nuz, phi_3d

    call make_grid(nxyz(3), .false., z, z_u)

    s_nuz = Solver_nuz(nxyz, real(Lxy, rp), real(z, rp), real(z_u, rp), bcs, &
                       PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)
    s_3d = Solver3D(nxyz, real([Lxy(1), Lxy(2), Lz], rp), bcs, &
                    PoisFFT_FiniteDifference2, planning=PoisFFT_Estimate)

    call random_number(tmp)
    rhs = real(tmp, rp)
    rhs = rhs - sum(rhs) / size(rhs)

    call Execute(s_nuz, phi_nuz, rhs)
    call Execute(s_3d, phi_3d, rhs)

    phi_nuz = phi_nuz - sum(phi_nuz) / size(phi_nuz)
    phi_3d = phi_3d - sum(phi_3d) / size(phi_3d)

    call check_err(bc_pair_name(bcs) // " FD2", rel_err(phi_nuz, phi_3d), tol)

    call Finalize(s_nuz)
    call Finalize(s_3d)
  end subroutine

end program test_nonuniform_z
