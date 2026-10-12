!Minimal assertion and reporting utilities for the PoisFFT unit tests.
!
!Every test program calls `ut_init` first and `ut_finish` last.
!`ut_finish` prints a summary and terminates with a non-zero exit code
!if any check failed, so that the test runner in SConstruct can detect it.
!
!Only failed checks are printed, unless the environment variable
!UT_VERBOSE is set to a non-empty value. Printing of failures is limited
!to UT_MAX_REPORT lines per program (default 40), the rest is only counted.
module unit_test_utils
  use iso_fortran_env, only: output_unit, real32, real64
  implicit none

  private
  public :: ut_init, ut_finish, check, check_err, rel_err, ut_section, &
            ut_failed_count

  integer, parameter :: sp = real32, dp = real64

  character(:), allocatable :: program_name
  character(:), allocatable :: section
  integer :: n_checks = 0, n_failed = 0, n_reported = 0
  integer :: max_report = 40
  logical :: verbose = .false.
  logical :: quiet = .false.

  interface rel_err
    module procedure rel_err_1d_sp, rel_err_2d_sp, rel_err_3d_sp
    module procedure rel_err_1d_dp, rel_err_2d_dp, rel_err_3d_dp
  end interface

  interface check_err
    module procedure check_err_sp, check_err_dp
  end interface

contains

  !silent=.true. suppresses all output, the checks are still counted
  !(used on the non-root MPI ranks, which evaluate the same global checks)
  subroutine ut_init(name, silent)
    character(*), intent(in) :: name
    logical, intent(in), optional :: silent
    character(64) :: val
    integer :: length, stat

    program_name = name
    section = ""
    n_checks = 0
    n_failed = 0
    n_reported = 0

    quiet = .false.
    if (present(silent)) quiet = silent

    call get_environment_variable("UT_VERBOSE", val, length, stat)
    verbose = (stat == 0 .and. length > 0) .and. .not. quiet

    call get_environment_variable("UT_MAX_REPORT", val, length, stat)
    if (stat == 0 .and. length > 0) read(val, *, iostat=stat) max_report
    if (quiet) max_report = 0

    if (.not. quiet) write(output_unit, '(a)') "=== " // name
  end subroutine


  subroutine ut_section(name)
    character(*), intent(in) :: name
    section = name
    if (verbose) write(output_unit, '(a)') "--- " // name
  end subroutine


  subroutine check(condition, label, detail)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    character(*), intent(in), optional :: detail
    character(:), allocatable :: line

    n_checks = n_checks + 1

    line = label
    if (len(section) > 0) line = section // ": " // line
    if (present(detail)) line = line // "  " // detail

    if (condition) then
      if (verbose) write(output_unit, '(a)') "  ok    " // line
    else
      n_failed = n_failed + 1
      n_reported = n_reported + 1
      if (quiet) then
        continue
      else if (n_reported <= max_report .or. verbose) then
        write(output_unit, '(a)') "  FAIL  " // line
      else if (n_reported == max_report + 1) then
        write(output_unit, '(a)') "  ...   further failures are only counted"
      end if
    end if
  end subroutine


  subroutine check_err_dp(label, err, tol)
    character(*), intent(in) :: label
    real(dp), intent(in) :: err, tol
    character(40) :: detail

    write(detail, '("err=",es9.2," tol=",es9.2)') err, tol
    ! a NaN error must fail as well, hence not (err > tol)
    call check(err <= tol, label, trim(detail))
  end subroutine

  subroutine check_err_sp(label, err, tol)
    character(*), intent(in) :: label
    real(sp), intent(in) :: err, tol

    call check_err_dp(label, real(err, dp), real(tol, dp))
  end subroutine


  integer function ut_failed_count() result(res)
    res = n_failed
  end function


  subroutine ut_finish()
    if (.not. quiet) &
      write(output_unit, '(a,i0,a,i0,a)') "=== " // program_name // ": ", &
            n_checks - n_failed, " of ", n_checks, " checks passed"
    flush(output_unit)
    if (n_failed > 0) error stop 1
  end subroutine


  !relative maximum-norm error max|a-b| / max|b|, the absolute error if b==0
  real(dp) function rel_err_1d_dp(a, b) result(res)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: scale
    scale = maxval(abs(b))
    if (scale == 0) scale = 1
    res = maxval(abs(a - b)) / scale
  end function

  real(dp) function rel_err_2d_dp(a, b) result(res)
    real(dp), intent(in) :: a(:,:), b(:,:)
    res = rel_err_1d_dp(reshape(a, [size(a)]), reshape(b, [size(b)]))
  end function

  real(dp) function rel_err_3d_dp(a, b) result(res)
    real(dp), intent(in) :: a(:,:,:), b(:,:,:)
    res = rel_err_1d_dp(reshape(a, [size(a)]), reshape(b, [size(b)]))
  end function

  real(sp) function rel_err_1d_sp(a, b) result(res)
    real(sp), intent(in) :: a(:), b(:)
    res = real(rel_err_1d_dp(real(a, dp), real(b, dp)), sp)
  end function

  real(sp) function rel_err_2d_sp(a, b) result(res)
    real(sp), intent(in) :: a(:,:), b(:,:)
    res = rel_err_1d_sp(reshape(a, [size(a)]), reshape(b, [size(b)]))
  end function

  real(sp) function rel_err_3d_sp(a, b) result(res)
    real(sp), intent(in) :: a(:,:,:), b(:,:,:)
    res = rel_err_1d_sp(reshape(a, [size(a)]), reshape(b, [size(b)]))
  end function

end module unit_test_utils
