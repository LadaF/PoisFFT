#include "precision.h"
!Unit test of the 2D solver, see test_eigen_1d.F90 for the method.
!The right-hand side is the tensor product of 1D eigenvectors in x and y,
!its eigenvalue is the sum of the 1D eigenvalues. All modes are checked
!for the combinations of boundary conditions in bc_combos.f90.
program test_eigen_2d
  use PoisFFT, Solver2D => SOLVER2D
  use reference_modes
  use unit_test_utils
  use bc_combos
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)
  real(dp), parameter :: Ls(2) = [2.3_dp, 1.7_dp]

  integer, parameter :: approxs(3) = [PoisFFT_Spectral, &
                                      PoisFFT_FiniteDifference2, &
                                      PoisFFT_FiniteDifference4]
  integer, parameter :: sizes(2,3) = reshape([6,7, 8,5, 3,4], [2,3])

  integer :: ic, ia, is

  call ut_init("test_eigen_2d " // PREC_NAME)

  do ic = 1, n_combos_2d
    do ia = 1, size(approxs)
      call ut_section("2D " // bc_pair_name(combos_2d(:,ic)) // &
                      " " // approx_name(approxs(ia)))
      do is = 1, size(sizes, 2)
        call test_all_modes(combos_2d(:,ic), approxs(ia), sizes(:,is))
      end do
    end do
  end do

  call ut_finish()

contains

  subroutine test_all_modes(bcs, approx, nxyz)
    integer, intent(in) :: bcs(4), approx, nxyz(2)
    type(Solver2D) :: solver
    real(rp) :: rhs(nxyz(1),nxyz(2)), phi(nxyz(1),nxyz(2)), exact(nxyz(1),nxyz(2))
    real(dp) :: vx(nxyz(1)), vy(nxyz(2)), lx, ly, lam, lam_min, err
    logical :: okx, oky
    integer :: mx, my, j
    character(60) :: label

    solver = Solver2D(nxyz, real(Ls, rp), bcs, approx, planning=PoisFFT_Estimate)
    lam_min = min_eigenvalue(bcs, nxyz, Ls, approx)

    do my = 1, nxyz(2)
      call mode_1d(bcs(3:4), nxyz(2), Ls(2), my, approx, vy, ly, oky)
      do mx = 1, nxyz(1)
        call mode_1d(bcs(1:2), nxyz(1), Ls(1), mx, approx, vx, lx, okx)
        lam = lx + ly
        write(label, '("n=",i0,"x",i0," mode ",i0,",",i0)') nxyz, mx, my

        do j = 1, nxyz(2)
          rhs(:,j) = real(vx * vy(j), rp)
        end do
        if (lam == 0) then
          exact = 0
        else
          do j = 1, nxyz(2)
            exact(:,j) = real(vx * vy(j) / lam, rp)
          end do
        end if

        call Execute(solver, phi, rhs)

        err = maxval(abs(phi - exact)) * lam_min / maxval(abs(rhs))
        call check_err(trim(label), err, real(tol, dp))
      end do
    end do

    call Finalize(solver)
  end subroutine

end program test_eigen_2d
