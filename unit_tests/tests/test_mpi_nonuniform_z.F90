#include "precision.h"
!Distributed version of test_nonuniform_z.F90 (compiled only with mpi=yes).
!
!  mpirun -n <npy*npz> test_mpi_nonuniform_z_dp <npy> <npz> [case]
!
!See test_nonuniform_z.F90 for the method and test_mpi_eigen_3d.F90 for the
!meaning of `case`. The z profile f is a fixed
!deterministic function, so that it is the same on all processes.
program test_mpi_nonuniform_z
  use mpi
  use PoisFFT, Solver_nuz => SOLVER3D_NUZ
  use reference_modes
  use unit_test_utils
  use mpi_test_grid
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 1000 * epsilon(1._rp)
  real(dp), parameter :: Lxy(2) = [2.3_dp, 1.7_dp], Lz = 1.1_dp
  integer, parameter :: gn(3) = [5, 6, 7]

  integer, parameter :: P = PoisFFT_Periodic, N = PoisFFT_Neumann, &
                        Ns = PoisFFT_NeumannStag, D = PoisFFT_Dirichlet, &
                        Ds = PoisFFT_DirichletStag
  !the combinations of test_nonuniform_z.F90
  integer, parameter :: n_combos = 12
  integer, parameter :: combos(6, n_combos) = reshape([ &
      P, P, P, P, Ns,Ns, &
      P, P, P, P, N, N,  &
      Ns,Ns,Ns,Ns,Ns,Ns, &
      N, N, N, N, N, N,  &
      P, P, D, D, N, N,  &
      D, D, P, P, Ns,Ns, &
      D, D, D, D, N, N,  &
      Ds,Ds,Ns,Ns,Ns,Ns, &
      N, D, P, P, N, N,  &
      Ns,Ds,Ds,Ns,Ns,Ns, &
      P, P, N, N, Ns,Ns, &
      N, N, Ns,Ns,N, N], [6, n_combos])
  integer, parameter :: approxs(2) = [PoisFFT_Spectral, PoisFFT_FiniteDifference2]
  integer :: ic, ia
  character(32) :: grid_name

  call grid_init(gn)

  write(grid_name, '(" grid ",i0,"x",i0)') npy, npz
  call ut_init("test_mpi_nonuniform_z " // PREC_NAME // trim(grid_name), silent=.not.is_root)

  do ic = 1, n_combos
    if (command_case /= 0 .and. command_case /= ic) cycle
    do ia = 1, size(approxs)
      call ut_section("MPI non-uniform z " // bc_pair_name(combos(:,ic)) // " " // &
                      approx_name(approxs(ia)))
      call test_modes(combos(:,ic), approxs(ia))
    end do
  end do

  call grid_finish()
  call ut_finish()

contains

  subroutine make_grid(nz, z, z_u)
    integer, intent(in) :: nz
    real(dp), intent(out) :: z(1:nz), z_u(0:nz)
    integer :: k

    do k = 0, nz
      z_u(k) = Lz * (real(k, dp) / nz)**1.5_dp
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


  real(dp) function global_mean(a) result(res)
    real(rp), intent(in) :: a(:,:,:)
    real(dp) :: s
    integer :: ie
    call MPI_Allreduce(real(sum(a), dp), s, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ie)
    res = s / product(gn)
  end function


  subroutine test_modes(bcs, approx)
    integer, intent(in) :: bcs(6), approx
    type(Solver_nuz) :: solver
    real(dp) :: z(gn(3)), z_u(0:gn(3)), f(gn(3)), af(gn(3))
    real(dp) :: vx(gn(1)), vy(gn(2)), lx, ly, err
    real(rp) :: scale
    real(rp), allocatable :: rhs(:,:,:), phi(:,:,:), exact(:,:,:)
    logical :: ok
    integer :: mx, my, i, j, k, kg
    character(80) :: label

    allocate(rhs(nxyz(1),nxyz(2),nxyz(3)), phi(nxyz(1),nxyz(2),nxyz(3)), &
             exact(nxyz(1),nxyz(2),nxyz(3)))

    call make_grid(gn(3), z, z_u)
    f = [(cos(1.3_dp * k) + 0.1_dp * k, k = 1, gn(3))]
    af = apply_z(f, z, z_u)

    !the solver gets the global z and z_u arrays
    solver = Solver_nuz(nxyz, real(Lxy, rp), real(z, rp), real(z_u, rp), bcs, approx, &
                        gnxyz=gn, offs=offs, mpi_comm=cart_comm, &
                        planning=PoisFFT_Estimate)

    do my = 1, gn(2)
      call mode_1d(bcs(3:4), gn(2), Lxy(2), my, approx, vy, ly, ok)
      do mx = 1, gn(1)
        call mode_1d(bcs(1:2), gn(1), Lxy(1), mx, approx, vx, lx, ok)

        do k = 1, nxyz(3)
          kg = k + offs(3)
          do j = 1, nxyz(2)
            do i = 1, nxyz(1)
              rhs(i,j,k) = real(vx(i+offs(1)) * vy(j+offs(2)) * (af(kg) + (lx + ly) * f(kg)), rp)
              exact(i,j,k) = real(vx(i+offs(1)) * vy(j+offs(2)) * f(kg), rp)
            end do
          end do
        end do

        call Execute(solver, phi, rhs)

        !see test_nonuniform_z.F90 for the scale
        scale = maxval(abs(exact))
        if (lx + ly == 0) then
          exact = exact - real(global_mean(exact), rp)
          phi = phi - real(global_mean(phi), rp)
        end if

        err = global_rel_err(real(maxval(abs(phi - exact)), dp), real(scale, dp))
        write(label, '("n=",i0,"x",i0,"x",i0," mode ",i0,",",i0)') gn, mx, my
        !see test_nonuniform_z.F90 for the tolerance of the x-y null mode
        call check_err(trim(label), err, real(merge(10*tol, tol, lx + ly == 0), dp))
      end do
    end do

    call Finalize(solver)
  end subroutine

end program test_mpi_nonuniform_z
