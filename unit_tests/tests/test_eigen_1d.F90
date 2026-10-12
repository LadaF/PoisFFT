#include "precision.h"
!Unit test of the 1D solver: eigenvalues, transform types and normalization.
!
!The right-hand side is a single discrete eigenvector v of the Laplacian with
!eigenvalue lam (from the independent module reference_modes), the exact
!solution is v / lam (0 for the null mode of periodic and Neumann problems).
!The error is normalized as described in reference_modes.f90.
!Every mode of every supported pair of boundary conditions is checked for the
!spectral, FD2 and FD4 approximations and for even and odd n.
program test_eigen_1d
  use PoisFFT, Solver1D => SOLVER1D
  use reference_modes
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)
  real(dp), parameter :: L = 2.3_dp

  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ds = PoisFFT_DirichletStag, &
                        Ns = PoisFFT_NeumannStag
  integer, parameter :: pairs(2,9) = reshape([P,P, D,D, Ds,Ds, N,N, Ns,Ns, &
                                              Ns,Ds, Ds,Ns, N,D, D,N], [2,9])
  integer, parameter :: approxs(3) = [PoisFFT_Spectral, &
                                      PoisFFT_FiniteDifference2, &
                                      PoisFFT_FiniteDifference4]
  integer, parameter :: sizes(5) = [2, 3, 8, 16, 17]

  integer :: ip, ia, is

  call ut_init("test_eigen_1d " // PREC_NAME)

  do ip = 1, size(pairs, 2)
    do ia = 1, size(approxs)
      call ut_section("1D " // bc_pair_name(pairs(:,ip)) // " " // approx_name(approxs(ia)))
      do is = 1, size(sizes)
        call test_all_modes(pairs(:,ip), approxs(ia), sizes(is))
      end do
    end do
  end do

  call ut_finish()

contains

  subroutine test_all_modes(bcs, approx, n)
    integer, intent(in) :: bcs(2), approx, n
    type(Solver1D) :: solver
    real(rp) :: rhs(n), phi(n), exact(n)
    real(dp) :: v(n), lam, lam_min, err
    logical :: ok
    integer :: m
    character(60) :: label

    solver = Solver1D([n], [real(L, rp)], bcs, approx, planning=PoisFFT_Estimate)
    lam_min = min_eigenvalue(bcs, [n], [L], approx)

    do m = 1, n
      call mode_1d(bcs, n, L, m, approx, v, lam, ok)
      write(label, '("n=",i0," mode ",i0)') n, m

      rhs = real(v, rp)
      if (lam == 0) then
        exact = 0
      else
        exact = real(v / lam, rp)
      end if

      call Execute(solver, phi, rhs)

      err = maxval(abs(phi - exact)) * lam_min / maxval(abs(rhs))
      call check_err(trim(label), err, real(tol, dp))
    end do

    call Finalize(solver)
  end subroutine

end program test_eigen_1d
