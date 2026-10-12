# PoisFFT unit tests

Self-contained test suite for the library in `../src`. The SConstruct in this
directory compiles the library sources itself (with runtime checks enabled), so
the tests do not depend on a previous build in `../bin` or `../lib`.

## Running

```
scons                      # build the serial library and tests, run all tests
scons mpi=yes              # MPI library, distributed tests on several process grids
scons build                # only build
scons tests=eigen_3d,fft   # run only the tests whose name contains one of the strings
scons prec=dp              # only double precision (default: both)
scons mpi=yes grids=2x2    # MPI process grids npy x npz (default 1x1,2x1,1x2,2x2)
scons src=/path/to/src     # test another source tree, e.g. a patched copy
```

Other options: `debug=no` (no `-fcheck`), `fpe=yes` (trap floating point
exceptions), `openmp=no`, `threads=N` (`OMP_NUM_THREADS` of the runs, default 2),
`timeout=S`, `fftw_inc=`, `fftw_lib=`, `pfft_inc=`, `pfft_lib=`. The compilers can
be set by the environment variables `FC` and `MPIFC`, the MPI launcher by `MPIRUN`.

The runner prints one line per test run (`PASS`, `FAIL`, `CRASH` or `TIMEOUT`)
and the failed checks; full logs are in `build/<variant>/logs`. Set `UT_VERBOSE=1`
to print every check, `UT_MAX_REPORT=N` to change the number of reported
failures per program (default 40). `scons` exits with an error if any run failed.

The MPI tests run all boundary-condition cases in one `mpirun` per process grid.
If such a run does not pass, the runner repeats it once per case, so that a crash
in one solver does not hide the results of the others (`split=yes` always runs the
cases separately).

## Layout

| File | Tests |
|---|---|
| `support/unit_test_utils.f90` | assertions (`check`, `check_err`), summary and exit code |
| `support/reference_modes.f90` | independent analytic eigenvectors and eigenvalues of the discrete 1D Laplacian for every boundary condition pair and approximation |
| `support/precision.h` | `-DSINGLE` selects the single-precision solver types |
| `support/bc_combos.f90` | the 39 combinations of boundary conditions of the 3D tests, shared by the serial and MPI tests |
| `support/mpi_test_grid.F90` | process grid, local block of the global grid, global error reduction |
| `tests/test_fft.F90` | `FFT_SP/DP`: all r2r kinds and the complex DFT in 1D/2D/3D against the FFTW definitions, order of the plan types, forward/backward normalization |
| `tests/test_eigen_1d.F90` | 1D solver: every mode, all 9 boundary-condition pairs, spectral/FD2/FD4, even and odd sizes |
| `tests/test_eigen_2d.F90` | 2D solver, same method |
| `tests/test_eigen_3d.F90` | 3D solver: every specialized solver routine (39 boundary-condition combinations) |
| `tests/test_solver_api.F90` | ghost layers, RHS not modified, repeated calls, `PoisFFT_Estimate` reproducibility, `PoisFFT_Measure`, linearity, zero mean of singular problems, re-creating solvers |
| `tests/test_nonuniform_z.F90` | `PoisFFT_Solver3D_nonuniform_z` against a manufactured solution on stretched and uniform grids, agreement with `PoisFFT_Solver3D` FD2 |
| `tests/test_error_reporting.F90` | errors returned in `ierr`/`errmsg` of the non-uniform z constructor: unsupported boundary conditions, too short or non-increasing `z`/`z_u`; finalizing such a solver |
| `tests/test_unsupported_bcs.F90` | without `ierr`, unsupported boundary conditions are rejected with "Unknown combination of boundary conditions." by all solver types, and an inconsistent z grid by the non-uniform z solver (serial, double precision only; every case runs as a separate process, the runner checks the output) |
| `tests/test_mpi_eigen_3d.F90` | distributed 3D solver, the same 39 combinations as the serial test |
| `tests/test_mpi_nonuniform_z.F90` | distributed non-uniform z solver |

## Method

Most tests use the right-hand side equal to one discrete eigenvector `v` of the
Laplacian (a tensor product of the 1D eigenvectors in 2D and 3D). The exact solution
is `v / lam`, or zero for the null mode of periodic and Neumann problems. Checking
every mode verifies the transform types, their normalization and every eigenvalue
independently. The error is normalized as `max|phi - exact| * lam_min / max|rhs|`,
which is `O(eps)` for a correct solver for all modes and grid sizes, so the
tolerance is `100 eps` in both precisions.

## Adding a test

Create `tests/test_<name>.F90` starting with `#include "precision.h"`, call
`ut_init`, the checks and `ut_finish`, and add it to `test_list` in `SConstruct`.
It is compiled in both precisions as `test_<name>_dp` and `test_<name>_sp`.
