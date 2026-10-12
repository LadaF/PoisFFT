!Independent reference for the PoisFFT unit tests.
!
!For every supported pair of boundary conditions in one direction it gives
!the discrete eigenvectors of the 1D Laplacian on the PoisFFT grid and the
!corresponding eigenvalues for the spectral, 2nd order and 4th order
!approximations. Nothing here uses the library code apart from the
!boundary condition and approximation constants.
!
!Grids (n unknowns, points j = 1..n, see usage_Fortran.txt):
!  Periodic        dx = L/n      x_j = (j-1) dx
!  Dirichlet       dx = L/(n+1)  x_j = j dx         (boundary points 0 and n+1)
!  Neumann         dx = L/(n-1)  x_j = (j-1) dx     (boundary points 1 and n)
!  DirichletStag,
!  NeumannStag     dx = L/n      x_j = (j-1/2) dx   (boundaries between points)
!  Neumann-Dirichlet  dx = L/n   x_j = (j-1) dx     (Dirichlet point n+1)
!  Dirichlet-Neumann  dx = L/n   x_j = j dx         (Dirichlet point 0)
!
!Mode m = 1..n. For a mode with the continuous wave number kappa and
!theta = kappa*dx/2 the eigenvalues are
!  spectral  -kappa**2
!  FD2       -(2 sin(theta) / dx)**2
!  FD4       -((27 sin(theta) - sin(3 theta)) / (12 dx))**2
module reference_modes
  use iso_fortran_env, only: dp => real64
  use PoisFFT_Parameters
  implicit none

  private
  public :: dp, pi, mode_1d, bc_supported_1d, bc_pair_name, approx_name, &
            ref_grid_dx, bc_letter, min_eigenvalue

  !The tests compare the solution of a mode with eigenvalue lam
  !using the error normalized as
  !   max|phi - exact| * lam_min / max|rhs|
  !where lam_min is the smallest non-zero |eigenvalue| of the problem.
  !Rounding errors of size eps*max|rhs| in the transforms of any mode
  !are amplified at most by 1/lam_min in the solution, so this error
  !is O(eps) independent of the mode and of the grid size.

  real(dp), parameter :: pi = 3.141592653589793238462643383279_dp

contains

  logical function bc_supported_1d(bcs) result(res)
    integer, intent(in) :: bcs(2)
    real(dp) :: v(1), lam
    logical :: ok
    call mode_1d(bcs, 1, 1._dp, 1, PoisFFT_FiniteDifference2, v, lam, ok, check_only=.true.)
    res = ok
  end function


  !The dx of the grid for the given boundary conditions.
  real(dp) function ref_grid_dx(bcs, n, L) result(dx)
    integer, intent(in) :: bcs(2), n
    real(dp), intent(in) :: L

    if (all(bcs == PoisFFT_Periodic)) then
      dx = L / n
    else if (all(bcs == PoisFFT_Dirichlet)) then
      dx = L / (n + 1)
    else if (all(bcs == PoisFFT_Neumann)) then
      dx = L / (n - 1)
    else
      dx = L / n
    end if
  end function


  !Eigenvector v(1:n) and eigenvalue lam of mode m (1..n).
  !ok is false for an unsupported pair of boundary conditions.
  subroutine mode_1d(bcs, n, L, m, approx, v, lam, ok, check_only)
    integer, intent(in) :: bcs(2), n, m, approx
    real(dp), intent(in) :: L
    real(dp), intent(out) :: v(:), lam
    logical, intent(out) :: ok
    logical, intent(in), optional :: check_only
    real(dp) :: dx, kappa, theta, x(n), h
    integer :: j, mm, kw
    logical :: only_check

    only_check = .false.
    if (present(check_only)) only_check = check_only

    ok = .true.
    v = 0
    lam = 0
    dx = ref_grid_dx(bcs, n, L)
    h = m - 0.5_dp

    if (all(bcs == PoisFFT_Periodic)) then
      if (only_check) return
      mm = m - 1
      kw = min(mm, n - mm)
      kappa = 2 * pi * kw / L
      x = [((j-1)*dx, j = 1, n)]
      if (mm <= n/2) then
        v = cos(2 * pi * mm * x / L)
      else
        v = sin(2 * pi * mm * x / L)
      end if

    else if (all(bcs == PoisFFT_Dirichlet)) then
      if (only_check) return
      kappa = pi * m / L
      x = [(j*dx, j = 1, n)]
      v = sin(kappa * x)

    else if (all(bcs == PoisFFT_DirichletStag)) then
      if (only_check) return
      kappa = pi * m / L
      x = [((j-0.5_dp)*dx, j = 1, n)]
      v = sin(kappa * x)

    else if (all(bcs == PoisFFT_Neumann)) then
      if (only_check) return
      kappa = pi * (m-1) / L
      x = [((j-1)*dx, j = 1, n)]
      v = cos(kappa * x)

    else if (all(bcs == PoisFFT_NeumannStag)) then
      if (only_check) return
      kappa = pi * (m-1) / L
      x = [((j-0.5_dp)*dx, j = 1, n)]
      v = cos(kappa * x)

    else if (all(bcs == [PoisFFT_NeumannStag, PoisFFT_DirichletStag])) then
      if (only_check) return
      kappa = pi * h / L
      x = [((j-0.5_dp)*dx, j = 1, n)]
      v = cos(kappa * x)

    else if (all(bcs == [PoisFFT_DirichletStag, PoisFFT_NeumannStag])) then
      if (only_check) return
      kappa = pi * h / L
      x = [((j-0.5_dp)*dx, j = 1, n)]
      v = sin(kappa * x)

    else if (all(bcs == [PoisFFT_Neumann, PoisFFT_Dirichlet])) then
      if (only_check) return
      kappa = pi * h / L
      x = [((j-1)*dx, j = 1, n)]
      v = cos(kappa * x)

    else if (all(bcs == [PoisFFT_Dirichlet, PoisFFT_Neumann])) then
      if (only_check) return
      kappa = pi * h / L
      x = [(j*dx, j = 1, n)]
      v = sin(kappa * x)

    else
      ok = .false.
      return
    end if

    theta = kappa * dx / 2

    select case (approx)
    case (PoisFFT_FiniteDifference2)
      lam = -(2 * sin(theta) / dx)**2
    case (PoisFFT_FiniteDifference4)
      lam = -((27 * sin(theta) - sin(3 * theta)) / (12 * dx))**2
    case default
      lam = -kappa**2
    end select
  end subroutine


  !The smallest non-zero |sum of 1D eigenvalues| over all modes of a 1D, 2D
  !or 3D problem, bcs(2*dims), nxyz(dims), Ls(dims).
  real(dp) function min_eigenvalue(bcs, nxyz, Ls, approx) result(res)
    integer, intent(in) :: bcs(:), nxyz(:), approx
    real(dp), intent(in) :: Ls(:)
    real(dp), allocatable :: lams(:,:), v(:)
    real(dp) :: lam
    integer :: d, i, j, k, nd, nmax
    integer :: ns(3)
    logical :: ok

    nd = size(nxyz)
    nmax = maxval(nxyz)
    allocate(lams(nmax, 3))
    lams = 0
    ns = 1
    ns(1:nd) = nxyz
    do d = 1, nd
      allocate(v(nxyz(d)))
      do i = 1, nxyz(d)
        call mode_1d(bcs(2*d-1:2*d), nxyz(d), Ls(d), i, approx, v, lams(i,d), ok)
      end do
      deallocate(v)
    end do

    res = huge(1._dp)
    do k = 1, ns(3)
      do j = 1, ns(2)
        do i = 1, ns(1)
          lam = abs(lams(i,1) + lams(j,2) + lams(k,3))
          if (lam > 0) res = min(res, lam)
        end do
      end do
    end do
  end function


  character(2) function bc_letter(bc) result(res)
    integer, intent(in) :: bc
    select case (bc)
    case (PoisFFT_Periodic);      res = "P"
    case (PoisFFT_Dirichlet);     res = "D"
    case (PoisFFT_Neumann);       res = "N"
    case (PoisFFT_DirichletStag); res = "Ds"
    case (PoisFFT_NeumannStag);   res = "Ns"
    case default;                 res = "?"
    end select
  end function


  function bc_pair_name(bcs) result(res)
    integer, intent(in) :: bcs(:)
    character(:), allocatable :: res
    integer :: i
    res = trim(bc_letter(bcs(1)))
    do i = 2, size(bcs)
      res = res // "," // trim(bc_letter(bcs(i)))
    end do
  end function


  function approx_name(approx) result(res)
    integer, intent(in) :: approx
    character(:), allocatable :: res
    select case (approx)
    case (PoisFFT_FiniteDifference2); res = "FD2"
    case (PoisFFT_FiniteDifference4); res = "FD4"
    case default;                     res = "spectral"
    end select
  end function

end module reference_modes
