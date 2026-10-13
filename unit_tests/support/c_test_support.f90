!C-callable access to the reference eigenmodes (reference_modes.f90) and to
!the reporting utilities (unit_test_utils.f90) for the tests of the C and
!C++ interface, declared in c_test_support.h. The strings are
!NUL-terminated C strings.
module c_test_support
  use iso_c_binding
  use reference_modes, only: dp, mode_1d, min_eigenvalue
  use unit_test_utils
  implicit none

  private

contains

  !Eigenvector v(0:n-1) and eigenvalue lam of mode m (1..n) for the pair of
  !boundary conditions bcs, returns 0 for an unsupported pair.
  integer(c_int) function ut_c_mode_1d(bcs, n, L, m, approx, v, lam) &
      result(res) bind(C, name="ut_c_mode_1d")
    integer(c_int), intent(in) :: bcs(2)
    integer(c_int), value :: n, m, approx
    real(c_double), value :: L
    real(c_double), intent(out) :: v(n), lam
    logical :: ok

    call mode_1d(int(bcs), int(n), real(L, dp), int(m), int(approx), v, lam, ok)
    res = merge(1, 0, ok)
  end function


  !The smallest non-zero |eigenvalue| of a problem with dims dimensions,
  !bcs(2*dims), nxyz(dims) and Ls(dims) in the same order of directions.
  real(c_double) function ut_c_min_eigenvalue(dims, bcs, nxyz, Ls, approx) &
      result(res) bind(C, name="ut_c_min_eigenvalue")
    integer(c_int), value :: dims, approx
    integer(c_int), intent(in) :: bcs(2*dims), nxyz(dims)
    real(c_double), intent(in) :: Ls(dims)

    res = min_eigenvalue(int(bcs), int(nxyz), real(Ls, dp), int(approx))
  end function


  subroutine ut_c_init(name) bind(C, name="ut_c_init")
    character(kind=c_char), intent(in) :: name(*)
    call ut_init(f_string(name))
  end subroutine


  subroutine ut_c_section(name) bind(C, name="ut_c_section")
    character(kind=c_char), intent(in) :: name(*)
    call ut_section(f_string(name))
  end subroutine


  subroutine ut_c_check(condition, label) bind(C, name="ut_c_check")
    integer(c_int), value :: condition
    character(kind=c_char), intent(in) :: label(*)
    call check(condition /= 0, f_string(label))
  end subroutine


  subroutine ut_c_check_err(label, err, tol) bind(C, name="ut_c_check_err")
    character(kind=c_char), intent(in) :: label(*)
    real(c_double), value :: err, tol
    call check_err(f_string(label), real(err, dp), real(tol, dp))
  end subroutine


  !terminates with a non-zero exit code if any check failed
  subroutine ut_c_finish() bind(C, name="ut_c_finish")
    call ut_finish()
  end subroutine


  function f_string(s) result(res)
    character(kind=c_char), intent(in) :: s(*)
    character(:), allocatable :: res
    integer :: n, i

    n = 0
    do while (s(n+1) /= c_null_char)
      n = n + 1
    end do
    allocate(character(n) :: res)
    do i = 1, n
      res(i:i) = s(i)
    end do
  end function

end module c_test_support
