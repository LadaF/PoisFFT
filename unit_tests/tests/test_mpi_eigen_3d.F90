#include "precision.h"
!Distributed version of test_eigen_3d.F90 (compiled only with mpi=yes).
!
!  mpirun -n <npy*npz> test_mpi_eigen_3d_dp <npy> <npz> [case]
!
!Every process builds its local block of the global eigenvector from the
!offsets of the PFFT decomposition, the error is reduced over all
!processes. `case` selects one combination of boundary conditions
!(1..n_combos), 0 or no argument runs all of them. The test runner starts
!each combination separately, so that a crash in one solver does not hide
!the results of the others.
program test_mpi_eigen_3d
  use PoisFFT, Solver3D => SOLVER3D
  use reference_modes
  use unit_test_utils
  use mpi_test_grid
  implicit none

  integer, parameter :: rp = RP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)
  real(dp), parameter :: Ls(3) = [2.3_dp, 1.7_dp, 1.1_dp]
  integer, parameter :: gn(3) = [5, 6, 7]

  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ds = PoisFFT_DirichletStag, &
                        Ns = PoisFFT_NeumannStag
  !one representative per specialized solver routine and its variants
  integer, parameter :: n_combos = 14
  integer, parameter :: combos(6, n_combos) = reshape([ &
      P, P, P, P, P, P,           & ! 1 FullPeriodic
      D, D, D, D, D, D,           & ! 2 FullDirichlet
      Ds,Ds,Ds,Ds,Ds,Ds,          & ! 3
      N, N, N, N, N, N,           & ! 4 FullNeumann
      Ns,Ns,Ns,Ns,Ns,Ns,          & ! 5
      P, P, P, P, Ns,Ns,          & ! 6 PPNs
      P, P, P, P, N, N,           & ! 7
      P, P, P, P, Ds,Ns,          & ! 8 PPD
      P, P, Ns,Ns,Ns,Ns,          & ! 9 PNsNs
      Ns,Ns,P, P, P, P,           & !10 NsPP
      P, P, Ns,Ns,P, P,           & !11 PNsP
      Ds,Ds,Ds,Ds,P, P,           & !12 DDP
      Ds,Ds,Ds,Ds,Ns,Ns,          & !13 2real1real
      Ns,Ns,Ds,Ds,Ds,Ds           & !14 1real2real
    ], [6, n_combos])
  integer, parameter :: approxs(3) = [PoisFFT_Spectral, &
                                      PoisFFT_FiniteDifference2, &
                                      PoisFFT_FiniteDifference4]
  integer :: ic, ia
  character(32) :: grid_name

  call grid_init(gn)

  write(grid_name, '(" grid ",i0,"x",i0)') npy, npz
  call ut_init("test_mpi_eigen_3d " // PREC_NAME // trim(grid_name), silent=.not.is_root)

  do ic = 1, n_combos
    if (command_case /= 0 .and. command_case /= ic) cycle
    do ia = 1, size(approxs)
      call ut_section("MPI 3D " // bc_pair_name(combos(:,ic)) // " " // approx_name(approxs(ia)))
      call test_all_modes(combos(:,ic), approxs(ia))
    end do
  end do

  call grid_finish()
  call ut_finish()

contains

  subroutine test_all_modes(bcs, approx)
    integer, intent(in) :: bcs(6), approx
    type(Solver3D) :: solver
    real(rp), allocatable :: rhs(:,:,:), phi(:,:,:), exact(:,:,:)
    real(dp) :: vx(gn(1)), vy(gn(2)), vz(gn(3)), lx, ly, lz, lam
    real(dp) :: err, lam_min
    logical :: ok(3)
    integer :: mx, my, mz, i, j, k
    character(80) :: label

    allocate(rhs(nxyz(1),nxyz(2),nxyz(3)), phi(nxyz(1),nxyz(2),nxyz(3)), &
             exact(nxyz(1),nxyz(2),nxyz(3)))

    solver = Solver3D(nxyz, real(Ls, rp), bcs, approx, &
                      gnxyz=gn, offs=offs, mpi_comm=cart_comm, &
                      planning=PoisFFT_Estimate)
    lam_min = min_eigenvalue(bcs, gn, Ls, approx)

    do mz = 1, gn(3)
      call mode_1d(bcs(5:6), gn(3), Ls(3), mz, approx, vz, lz, ok(3))
      do my = 1, gn(2)
        call mode_1d(bcs(3:4), gn(2), Ls(2), my, approx, vy, ly, ok(2))
        do mx = 1, gn(1)
          call mode_1d(bcs(1:2), gn(1), Ls(1), mx, approx, vx, lx, ok(1))
          lam = lx + ly + lz

          do k = 1, nxyz(3)
            do j = 1, nxyz(2)
              do i = 1, nxyz(1)
                rhs(i,j,k) = real(vx(i+offs(1)) * vy(j+offs(2)) * vz(k+offs(3)), rp)
                if (lam == 0) then
                  exact(i,j,k) = 0
                else
                  exact(i,j,k) = real(vx(i+offs(1)) * vy(j+offs(2)) * vz(k+offs(3)) / lam, rp)
                end if
              end do
            end do
          end do

          call Execute(solver, phi, rhs)

          !normalized as described in reference_modes.f90
          err = global_rel_err(real(maxval(abs(phi - exact)), dp) * lam_min, &
                               real(maxval(abs(rhs)), dp))
          write(label, '("n=",i0,"x",i0,"x",i0," mode ",i0,",",i0,",",i0)') gn, mx, my, mz
          call check_err(trim(label), err, real(tol, dp))
        end do
      end do
    end do

    call Finalize(solver)
  end subroutine

end program test_mpi_eigen_3d
