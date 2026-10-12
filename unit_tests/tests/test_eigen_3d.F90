#include "precision.h"
!Unit test of the 3D solver, see test_eigen_1d.F90 for the method.
!The right-hand side is the tensor product of 1D eigenvectors in x, y and z,
!its eigenvalue is the sum of the 1D eigenvalues. All modes are checked for
!every combination of boundary conditions for which PoisFFT_Solver3D has a
!specialized solver, so that each solver routine and its splitting into
!1D and 2D transforms is exercised.
program test_eigen_3d
  use PoisFFT, Solver3D => SOLVER3D
  use reference_modes
  use unit_test_utils
  use bc_combos
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)
  real(dp), parameter :: Ls(3) = [2.3_dp, 1.7_dp, 1.1_dp]

  integer, parameter :: approxs(3) = [PoisFFT_Spectral, &
                                      PoisFFT_FiniteDifference2, &
                                      PoisFFT_FiniteDifference4]
  integer, parameter :: sizes(3,2) = reshape([4,5,6, 6,4,5], [3,2])

  integer :: ic, ia, is

  call ut_init("test_eigen_3d " // PREC_NAME)

  do ic = 1, n_combos_3d
    do ia = 1, size(approxs)
      call ut_section("3D " // bc_pair_name(combos_3d(:,ic)) // " " // approx_name(approxs(ia)))
      do is = 1, size(sizes, 2)
        call test_all_modes(combos_3d(:,ic), approxs(ia), sizes(:,is))
      end do
    end do
  end do

  call ut_finish()

contains

  subroutine test_all_modes(bcs, approx, nxyz)
    integer, intent(in) :: bcs(6), approx, nxyz(3)
    type(Solver3D) :: solver
    real(rp), allocatable :: rhs(:,:,:), phi(:,:,:), exact(:,:,:)
    real(dp) :: vx(nxyz(1)), vy(nxyz(2)), vz(nxyz(3)), lx, ly, lz, lam, lam_min, err
    logical :: ok(3)
    integer :: mx, my, mz, j, k
    character(80) :: label

    allocate(rhs(nxyz(1),nxyz(2),nxyz(3)), phi(nxyz(1),nxyz(2),nxyz(3)), exact(nxyz(1),nxyz(2),nxyz(3)))

    solver = Solver3D(nxyz, real(Ls, rp), bcs, approx, planning=PoisFFT_Estimate)
    lam_min = min_eigenvalue(bcs, nxyz, Ls, approx)

    do mz = 1, nxyz(3)
      call mode_1d(bcs(5:6), nxyz(3), Ls(3), mz, approx, vz, lz, ok(3))
      do my = 1, nxyz(2)
        call mode_1d(bcs(3:4), nxyz(2), Ls(2), my, approx, vy, ly, ok(2))
        do mx = 1, nxyz(1)
          call mode_1d(bcs(1:2), nxyz(1), Ls(1), mx, approx, vx, lx, ok(1))
          lam = lx + ly + lz
          write(label, '("n=",i0,"x",i0,"x",i0," mode ",i0,",",i0,",",i0)') nxyz, mx, my, mz

          do k = 1, nxyz(3)
            do j = 1, nxyz(2)
              rhs(:,j,k) = real(vx * vy(j) * vz(k), rp)
              if (lam == 0) then
                exact(:,j,k) = 0
              else
                exact(:,j,k) = real(vx * vy(j) * vz(k) / lam, rp)
              end if
            end do
          end do

          call Execute(solver, phi, rhs)

          err = maxval(abs(phi - exact)) * lam_min / maxval(abs(rhs))
          call check_err(trim(label), err, real(tol, dp))
        end do
      end do
    end do

    call Finalize(solver)
  end subroutine

end program test_eigen_3d
