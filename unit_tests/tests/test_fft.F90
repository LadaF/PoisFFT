#include "precision.h"
!Unit test of the FFT wrapper modules FFT_SP / FFT_DP (fft-inc.f90).
!
!Every real-to-real transform kind and the complex DFT is created through
!PoisFFT_Plan1D/2D/3D and compared with a direct O(n**2) evaluation of the
!definition from the FFTW manual. The multidimensional checks also verify
!the order in which the transform kinds are applied: plantypes(1) belongs
!to the slowest varying (last Fortran) dimension, because the plans are
!created with the dimensions in the C order.
!Finally, the forward/backward pairs used by the solvers must give the
!identity multiplied by the solver normalization factor.
program test_fft
  use FFT_MODULE
  use PoisFFT_Parameters, only: PoisFFT_Estimate
  use reference_modes, only: dp, pi
  use unit_test_utils
  implicit none

  integer, parameter :: rp = RP_KIND, cp = CP_KIND
  real(rp), parameter :: tol = 100 * epsilon(1._rp)

  integer, parameter :: kinds(8) = [FFT_RealEven00, FFT_RealEven01, &
                                    FFT_RealEven10, FFT_RealEven11, &
                                    FFT_RealOdd00,  FFT_RealOdd01, &
                                    FFT_RealOdd10,  FFT_RealOdd11]
  character(8), parameter :: kind_names(8) = [character(8) :: "REDFT00", "REDFT01", &
                                    "REDFT10", "REDFT11", "RODFT00", "RODFT01", &
                                    "RODFT10", "RODFT11"]
  integer :: ik, n, i

  call ut_init("test_fft " // PREC_NAME)

  call ut_section("1D real-to-real")
  do ik = 1, size(kinds)
    do n = 2, 9
      call test_r2r_1d(ik, n)
    end do
  end do

  call ut_section("1D complex")
  do n = 1, 9
    call test_c2c_1d(n, FFTW_FORWARD)
    call test_c2c_1d(n, FFTW_BACKWARD)
  end do

  call ut_section("2D real-to-real, plantypes(1) acts on dimension 2")
  do i = 1, size(kinds)
    call test_r2r_2d(i, 1 + mod(i+2, size(kinds)), [5, 6])
  end do

  call ut_section("2D complex")
  call test_c2c_2d([5, 6], FFTW_FORWARD)
  call test_c2c_2d([4, 3], FFTW_BACKWARD)

  call ut_section("3D real-to-real, plantypes(1) acts on dimension 3")
  call test_r2r_3d([1, 6, 3], [4, 5, 3])
  call test_r2r_3d([8, 2, 7], [3, 4, 5])

  call ut_section("3D complex")
  call test_c2c_3d([3, 4, 5], FFTW_FORWARD)
  call test_c2c_3d([4, 3, 2], FFTW_BACKWARD)

  call ut_section("forward/backward pairs of the solvers")
  do n = 2, 9
    call test_round_trip(FFT_RealOdd00,  FFT_RealOdd00,  n, 2*(n+1), "Dirichlet")
    call test_round_trip(FFT_RealOdd10,  FFT_RealOdd01,  n, 2*n,     "DirichletStag")
    call test_round_trip(FFT_RealEven00, FFT_RealEven00, n, 2*(n-1), "Neumann")
    call test_round_trip(FFT_RealEven10, FFT_RealEven01, n, 2*n,     "NeumannStag")
    call test_round_trip(FFT_RealEven11, FFT_RealEven11, n, 2*n,     "NeumannStag-DirichletStag")
    call test_round_trip(FFT_RealOdd11,  FFT_RealOdd11,  n, 2*n,     "DirichletStag-NeumannStag")
    call test_round_trip(FFT_RealEven01, FFT_RealEven10, n, 2*n,     "Neumann-Dirichlet")
    call test_round_trip(FFT_RealOdd01,  FFT_RealOdd10,  n, 2*n,     "Dirichlet-Neumann")
  end do

  call ut_finish()

contains

  !Direct evaluation of the FFTW real-to-real transforms (FFTW manual, 4.8.3).
  function naive_r2r(kind, x) result(y)
    integer, intent(in) :: kind
    real(dp), intent(in) :: x(0:)
    real(dp) :: y(0:size(x)-1)
    integer :: j, k, n

    n = size(x)
    do k = 0, n-1
      select case (kind)
      case (FFT_RealEven00)
        y(k) = x(0) + (-1)**k * x(n-1)
        do j = 1, n-2
          y(k) = y(k) + 2 * x(j) * cos(pi * j * k / (n-1))
        end do
      case (FFT_RealEven10)
        y(k) = 2 * sum([(x(j) * cos(pi * (j+0.5_dp) * k / n), j = 0, n-1)])
      case (FFT_RealEven01)
        y(k) = x(0) + 2 * sum([(x(j) * cos(pi * j * (k+0.5_dp) / n), j = 1, n-1)])
      case (FFT_RealEven11)
        y(k) = 2 * sum([(x(j) * cos(pi * (j+0.5_dp) * (k+0.5_dp) / n), j = 0, n-1)])
      case (FFT_RealOdd00)
        y(k) = 2 * sum([(x(j) * sin(pi * (j+1) * (k+1) / (n+1._dp)), j = 0, n-1)])
      case (FFT_RealOdd10)
        y(k) = 2 * sum([(x(j) * sin(pi * (j+0.5_dp) * (k+1) / n), j = 0, n-1)])
      case (FFT_RealOdd01)
        y(k) = (-1)**k * x(n-1) + &
               2 * sum([(x(j) * sin(pi * (j+1) * (k+0.5_dp) / n), j = 0, n-2)])
      case (FFT_RealOdd11)
        y(k) = 2 * sum([(x(j) * sin(pi * (j+0.5_dp) * (k+0.5_dp) / n), j = 0, n-1)])
      end select
    end do
  end function


  function naive_dft(x, dir) result(y)
    complex(dp), intent(in) :: x(0:)
    integer, intent(in) :: dir
    complex(dp) :: y(0:size(x)-1)
    integer :: j, k, n
    real(dp) :: s

    n = size(x)
    s = merge(-1._dp, 1._dp, dir == FFTW_FORWARD)
    do k = 0, n-1
      y(k) = sum([(x(j) * exp(cmplx(0, s * 2 * pi * j * k / n, dp)), j = 0, n-1)])
    end do
  end function


  subroutine fill_random(x)
    real(dp), intent(out) :: x(:)
    call random_number(x)
    x = x - 0.5_dp
  end subroutine


  subroutine test_r2r_1d(ik, n)
    integer, intent(in) :: ik, n
    type(PoisFFT_Solver1D) :: D
    type(PoisFFT_Plan1D) :: plan
    real(dp) :: x(n), y(n)
    character(40) :: label

    D%nx = n; D%gnx = n; D%nxyz = [n]; D%cnt = n; D%gcnt = n
    D%planning = PoisFFT_Estimate
    call allocate_fftw_real(D)
    plan = PoisFFT_Plan1D(D, [kinds(ik)])

    call fill_random(x)
    D%rwork = real(x, rp)
    call Execute(plan, D%rwork)
    y = naive_r2r(kinds(ik), real(real(x, rp), dp))

    write(label, '(a," n=",i0)') trim(kind_names(ik)), n
    call check_err(trim(label), rel_err(real(D%rwork, dp), y), real(tol, dp) * n)

    call Finalize(plan)
    call deallocate_fftw(D%rwork)
  end subroutine


  subroutine test_c2c_1d(n, dir)
    integer, intent(in) :: n, dir
    type(PoisFFT_Solver1D) :: D
    type(PoisFFT_Plan1D) :: plan
    real(dp) :: re(n), im(n)
    complex(dp) :: x(n), y(n)
    character(40) :: label

    D%nx = n; D%gnx = n; D%nxyz = [n]; D%cnt = n; D%gcnt = n
    D%planning = PoisFFT_Estimate
    call allocate_fftw_complex(D)
    plan = PoisFFT_Plan1D(D, [FFT_Complex, dir])

    call fill_random(re)
    call fill_random(im)
    x = cmplx(re, im, dp)
    D%cwork = cmplx(x, kind=cp)
    y = naive_dft(cmplx(D%cwork, kind=dp), dir)
    call Execute(plan, D%cwork)

    write(label, '("DFT dir=",i0," n=",i0)') dir, n
    call check_err(trim(label), &
                   max(rel_err(real(D%cwork, dp), real(y)), rel_err(real(aimag(D%cwork), dp), aimag(y))), &
                   real(tol, dp) * max(n, 1))

    call Finalize(plan)
    call deallocate_fftw(D%cwork)
  end subroutine


  !Kinds k1 = plantypes(1) and k2 = plantypes(2). Expected: k2 along
  !dimension 1 (x, fastest) and k1 along dimension 2 (y).
  subroutine test_r2r_2d(ik1, ik2, ns)
    integer, intent(in) :: ik1, ik2, ns(2)
    type(PoisFFT_Solver2D) :: D
    type(PoisFFT_Plan2D) :: plan
    real(dp) :: x(ns(1),ns(2)), y(ns(1),ns(2))
    integer :: i, j
    character(60) :: label

    D%nx = ns(1); D%ny = ns(2); D%gnx = ns(1); D%gny = ns(2)
    D%nxyz = ns; D%cnt = product(ns); D%gcnt = product(ns)
    D%planning = PoisFFT_Estimate
    call allocate_fftw_real(D)
    plan = PoisFFT_Plan2D(D, [kinds(ik1), kinds(ik2)])

    do j = 1, ns(2)
      call fill_random(x(:,j))
    end do
    x = real(real(x, rp), dp)
    D%rwork = real(x, rp)
    call Execute(plan, D%rwork)

    y = x
    do j = 1, ns(2)
      y(:,j) = naive_r2r(kinds(ik2), y(:,j))
    end do
    do i = 1, ns(1)
      y(i,:) = naive_r2r(kinds(ik1), y(i,:))
    end do

    write(label, '("[",a,",",a,"]")') trim(kind_names(ik1)), trim(kind_names(ik2))
    call check_err(trim(label), rel_err(real(D%rwork, dp), y), real(tol, dp) * product(ns))

    call Finalize(plan)
    call deallocate_fftw(D%rwork)
  end subroutine


  subroutine test_c2c_2d(ns, dir)
    integer, intent(in) :: ns(2), dir
    type(PoisFFT_Solver2D) :: D
    type(PoisFFT_Plan2D) :: plan
    real(dp) :: re(ns(1),ns(2)), im(ns(1),ns(2))
    complex(dp) :: y(ns(1),ns(2))
    integer :: i, j
    character(40) :: label

    D%nx = ns(1); D%ny = ns(2); D%gnx = ns(1); D%gny = ns(2)
    D%nxyz = ns; D%cnt = product(ns); D%gcnt = product(ns)
    D%planning = PoisFFT_Estimate
    call allocate_fftw_complex(D)
    plan = PoisFFT_Plan2D(D, [FFT_Complex, dir])

    do j = 1, ns(2)
      call fill_random(re(:,j))
      call fill_random(im(:,j))
    end do
    D%cwork = cmplx(re, im, cp)
    y = cmplx(D%cwork, kind=dp)
    call Execute(plan, D%cwork)

    do j = 1, ns(2)
      y(:,j) = naive_dft(y(:,j), dir)
    end do
    do i = 1, ns(1)
      y(i,:) = naive_dft(y(i,:), dir)
    end do

    write(label, '("DFT dir=",i0," n=",i0,"x",i0)') dir, ns
    call check_err(trim(label), &
                   max(rel_err(real(D%cwork, dp), real(y)), rel_err(real(aimag(D%cwork), dp), aimag(y))), &
                   real(tol, dp) * product(ns))

    call Finalize(plan)
    call deallocate_fftw(D%cwork)
  end subroutine


  !iks = plantypes(1:3); expected kinds(iks(3)) along x, iks(2) along y,
  !iks(1) along z
  subroutine test_r2r_3d(iks, ns)
    integer, intent(in) :: iks(3), ns(3)
    type(PoisFFT_Solver3D) :: D
    type(PoisFFT_Plan3D) :: plan
    real(dp) :: x(ns(1),ns(2),ns(3)), y(ns(1),ns(2),ns(3))
    integer :: i, j, k
    character(60) :: label

    D%nx = ns(1); D%ny = ns(2); D%nz = ns(3)
    D%gnx = ns(1); D%gny = ns(2); D%gnz = ns(3)
    D%nxyz = ns; D%cnt = product(ns); D%gcnt = product(ns)
    D%planning = PoisFFT_Estimate
    call allocate_fftw_real(D)
    plan = PoisFFT_Plan3D(D, kinds(iks))

    do k = 1, ns(3)
      do j = 1, ns(2)
        call fill_random(x(:,j,k))
      end do
    end do
    x = real(real(x, rp), dp)
    D%rwork = real(x, rp)
    call Execute(plan, D%rwork)

    y = x
    do k = 1, ns(3)
      do j = 1, ns(2)
        y(:,j,k) = naive_r2r(kinds(iks(3)), y(:,j,k))
      end do
    end do
    do k = 1, ns(3)
      do i = 1, ns(1)
        y(i,:,k) = naive_r2r(kinds(iks(2)), y(i,:,k))
      end do
    end do
    do j = 1, ns(2)
      do i = 1, ns(1)
        y(i,j,:) = naive_r2r(kinds(iks(1)), y(i,j,:))
      end do
    end do

    write(label, '("[",a,",",a,",",a,"]")') (trim(kind_names(iks(i))), i = 1, 3)
    call check_err(trim(label), rel_err(real(D%rwork, dp), y), real(tol, dp) * product(ns))

    call Finalize(plan)
    call deallocate_fftw(D%rwork)
  end subroutine


  subroutine test_c2c_3d(ns, dir)
    integer, intent(in) :: ns(3), dir
    type(PoisFFT_Solver3D) :: D
    type(PoisFFT_Plan3D) :: plan
    real(dp) :: re(ns(1),ns(2),ns(3)), im(ns(1),ns(2),ns(3))
    complex(dp) :: y(ns(1),ns(2),ns(3))
    integer :: i, j, k
    character(40) :: label

    D%nx = ns(1); D%ny = ns(2); D%nz = ns(3)
    D%gnx = ns(1); D%gny = ns(2); D%gnz = ns(3)
    D%nxyz = ns; D%cnt = product(ns); D%gcnt = product(ns)
    D%planning = PoisFFT_Estimate
    call allocate_fftw_complex(D)
    plan = PoisFFT_Plan3D(D, [FFT_Complex, dir])

    do k = 1, ns(3)
      do j = 1, ns(2)
        call fill_random(re(:,j,k))
        call fill_random(im(:,j,k))
      end do
    end do
    D%cwork = cmplx(re, im, cp)
    y = cmplx(D%cwork, kind=dp)
    call Execute(plan, D%cwork)

    do k = 1, ns(3)
      do j = 1, ns(2)
        y(:,j,k) = naive_dft(y(:,j,k), dir)
      end do
    end do
    do k = 1, ns(3)
      do i = 1, ns(1)
        y(i,:,k) = naive_dft(y(i,:,k), dir)
      end do
    end do
    do j = 1, ns(2)
      do i = 1, ns(1)
        y(i,j,:) = naive_dft(y(i,j,:), dir)
      end do
    end do

    write(label, '("DFT dir=",i0," n=",i0,"x",i0,"x",i0)') dir, ns
    call check_err(trim(label), &
                   max(rel_err(real(D%cwork, dp), real(y)), rel_err(real(aimag(D%cwork), dp), aimag(y))), &
                   real(tol, dp) * product(ns))

    call Finalize(plan)
    call deallocate_fftw(D%cwork)
  end subroutine


  subroutine test_round_trip(kf, kb, n, norm, name)
    integer, intent(in) :: kf, kb, n, norm
    character(*), intent(in) :: name
    type(PoisFFT_Solver1D) :: D
    type(PoisFFT_Plan1D) :: forward, backward
    real(dp) :: x(n)
    character(60) :: label

    D%nx = n; D%gnx = n; D%nxyz = [n]; D%cnt = n; D%gcnt = n
    D%planning = PoisFFT_Estimate
    call allocate_fftw_real(D)
    forward = PoisFFT_Plan1D(D, [kf])
    backward = PoisFFT_Plan1D(D, [kb])

    call fill_random(x)
    x = real(real(x, rp), dp)
    D%rwork = real(x, rp)
    call Execute(forward, D%rwork)
    call Execute(backward, D%rwork)

    write(label, '(a," n=",i0)') name, n
    call check_err(trim(label), rel_err(real(D%rwork, dp) / norm, x), real(tol, dp) * n)

    call Finalize(forward)
    call Finalize(backward)
    call deallocate_fftw(D%rwork)
  end subroutine

end program test_fft
