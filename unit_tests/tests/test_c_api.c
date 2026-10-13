/* Unit test of the C interface (poisfft.h, c_binding.c, c_binding.f90).

     test_c_api_dp                              serial build
     mpirun -n <npy*npz> test_mpi_c_api_dp <npy> <npz> [case]   MPI build

   The solvers themselves are tested from Fortran, this test checks the
   translation between C and Fortran: the C row-major order of the arrays
   and of the arguments (sizes, lengths and boundary conditions, a different
   size and pair in every direction), the ghost points (different numbers
   in every direction and for Phi and RHS), the approximation and the number
   of threads, for 1D, 2D and 3D. Compiled with -DSINGLE for the float
   functions.

   With MPI it also checks the global sizes, the offsets and the
   communicator. The 3D solvers use the 2D PFFT process grid npy x npz of
   the command line, the 1D and 2D solvers a 1D communicator of all
   processes. The local blocks of the 2D and 3D solvers come from PFFT, the
   1D grid is split into contiguous blocks. The case argument is not used.

   The right-hand side is a product of 1D eigenvectors from the Fortran
   reference (c_test_support.f90), the error is normalized as in the
   Fortran tests. The ghost points of Phi must not be changed by the
   solver, those of RHS contain a value that would spoil the solution. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#ifdef MPI
#include <mpi.h>
#include <pfft.h>
#endif
#include "poisfft.h"
#include "c_test_support.h"

#ifdef SINGLE
typedef float real;
typedef poisfft_solver_f solver_t;
#define solver_new poisfft_solver_f_new
#define solver_execute poisfft_solver_f_execute
#define solver_finalize poisfft_solver_f_finalize
#define REAL_EPS FLT_EPSILON
#define PREC_NAME "SP"
#else
typedef double real;
typedef poisfft_solver solver_t;
#define solver_new poisfft_solver_new
#define solver_execute poisfft_solver_execute
#define solver_finalize poisfft_solver_finalize
#define REAL_EPS DBL_EPSILON
#define PREC_NAME "DP"
#endif

static const double Ls[3] = {2.3, 1.7, 1.1};
/* the ghost points of RHS and Phi */
static const real ghost_rhs = -777, ghost_phi = 12345;
/* the numbers of ghost points of Phi and RHS in the directions 0, 1, 2 */
static const int ng_phi[3] = {1, 0, 2}, ng_rhs[3] = {0, 2, 1};

#ifdef MPI
/* the communicators of the 1D and 2D solvers and of the 3D solvers */
static MPI_Comm comm_2d, comm_3d;
#endif

static const char *bc_name(int bc)
{
  if (bc == POISFFT_PERIODIC) return "P";
  if (bc == POISFFT_DIRICHLET) return "D";
  if (bc == POISFFT_NEUMANN) return "N";
  if (bc == POISFFT_DIRICHLET_STAG) return "Ds";
  if (bc == POISFFT_NEUMANN_STAG) return "Ns";
  return "?";
}

static const char *approx_name(int approx)
{
  if (approx == POISFFT_FINITE_DIFFERENCE_2) return "FD2";
  if (approx == POISFFT_FINITE_DIFFERENCE_4) return "FD4";
  return "spectral";
}

/* maximum over all processes, a NaN gives infinity */
static double global_max(double x)
{
  if (x != x) x = HUGE_VAL;
#ifdef MPI
  MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return x;
}

/* true on all processes */
static int global_all(int x)
{
#ifdef MPI
  MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
  return x;
}

/* index of the point (i[0], i[1], i[2]) of an array with ng ghost points
   and sizes s, i is the local index of the unknown (0..n-1) */
static size_t idx(const int *i, const int *ng, const int *s)
{
  return ((size_t)(i[0] + ng[0]) * s[1] + (i[1] + ng[1])) * s[2] + (i[2] + ng[2]);
}

/* All modes of one problem. dims directions with the global sizes gns,
   boundary conditions bcs (2*dims) in the C order. Internally the
   directions are padded to three at the front, the slowest varying index,
   with n=1. */
static void test_problem(int dims, const int *gns, const int *bcs, int approx,
                         int ghosts)
{
  int ns[3], offs[3];
  int n[3], gn[3], off[3], gp[3], gr[3], sp[3], sr[3], m[3], i[3], d, pd;
  real Lr[3];
  double *v[3], l[3], lam, lam_min, err, rhs_max;
  real *phi, *rhs;
  size_t np, nr, k;
  int ghost_ok;
  char section[120], label[120];
  solver_t S;

  /* the local block */
#ifdef MPI
  if (dims == 1) {
    int rank, nprocs;
    MPI_Comm_rank(comm_2d, &rank);
    MPI_Comm_size(comm_2d, &nprocs);
    ns[0] = gns[0] / nprocs + (rank < gns[0] % nprocs);
    offs[0] = rank * (gns[0] / nprocs) + (rank < gns[0] % nprocs ? rank : gns[0] % nprocs);
  } else {
    MPI_Comm comm = dims == 3 ? comm_3d : comm_2d;
    ptrdiff_t gn_p[3], ln[3], lo[3], ln2[3], lo2[3];
    for (d = 0; d < dims; d++) gn_p[d] = gns[d];
    pfft_local_size_dft(dims, gn_p, comm, PFFT_TRANSPOSED_NONE, ln, lo, ln2, lo2);
    for (d = 0; d < dims; d++) {
      ns[d] = (int) ln[d];
      offs[d] = (int) lo[d];
    }
  }
#else
  for (d = 0; d < dims; d++) {
    ns[d] = gns[d];
    offs[d] = 0;
  }
#endif

  for (d = 0; d < 3; d++) {
    n[d] = 1; gn[d] = 1; off[d] = 0; gp[d] = 0; gr[d] = 0;
  }
  for (d = 0; d < dims; d++) {
    pd = 3 - dims + d;
    n[pd] = ns[d];
    gn[pd] = gns[d];
    off[pd] = offs[d];
    if (ghosts) {
      gp[pd] = ng_phi[d];
      gr[pd] = ng_rhs[d];
    }
    Lr[d] = (real) Ls[d];
  }
  for (d = 0; d < 3; d++) {
    sp[d] = n[d] + 2 * gp[d];
    sr[d] = n[d] + 2 * gr[d];
    v[d] = malloc(gn[d] * sizeof(double));
    v[d][0] = 1;
    l[d] = 0;
  }
  np = (size_t) sp[0] * sp[1] * sp[2];
  nr = (size_t) sr[0] * sr[1] * sr[2];
  phi = malloc(np * sizeof(real));
  rhs = malloc(nr * sizeof(real));

  k = sprintf(section, "C %dD", dims);
  for (d = 0; d < dims; d++)
    k += sprintf(section + k, " %s,%s", bc_name(bcs[2*d]), bc_name(bcs[2*d+1]));
  sprintf(section + k, " %s%s", approx_name(approx), ghosts ? " ghosts" : "");
  ut_c_section(section);

#ifdef MPI
  S = solver_new(dims, ns, Lr, bcs, approx, gns, offs,
                 dims == 3 ? &comm_3d : &comm_2d, ghosts ? 2 : 0);
#else
  S = solver_new(dims, ns, Lr, bcs, approx, NULL, NULL, NULL, ghosts ? 2 : 0);
#endif
  ut_c_check(global_all(S.D != NULL), "solver created");
  if (S.D == NULL) return;

  lam_min = ut_c_min_eigenvalue(dims, bcs, gns, Ls, approx);

  for (m[0] = 1; m[0] <= gn[0]; m[0]++)
  for (m[1] = 1; m[1] <= gn[1]; m[1]++)
  for (m[2] = 1; m[2] <= gn[2]; m[2]++) {
    for (d = 0; d < dims; d++) {
      pd = 3 - dims + d;
      ut_c_mode_1d(&bcs[2*d], gns[d], Ls[d], m[pd], approx, v[pd], &l[pd]);
    }
    lam = l[0] + l[1] + l[2];

    for (k = 0; k < nr; k++) rhs[k] = ghost_rhs;
    for (k = 0; k < np; k++) phi[k] = ghost_phi;
    for (i[0] = 0; i[0] < n[0]; i[0]++)
    for (i[1] = 0; i[1] < n[1]; i[1]++)
    for (i[2] = 0; i[2] < n[2]; i[2]++)
      rhs[idx(i, gr, sr)] = (real) (v[0][i[0]+off[0]] * v[1][i[1]+off[1]] * v[2][i[2]+off[2]]);

    if (ghosts) {
      /* only the first dims numbers, in the C order of the directions */
      solver_execute(S, phi, rhs, ng_phi, ng_rhs);
    } else {
      solver_execute(S, phi, rhs, NULL, NULL);
    }

    err = 0;
    rhs_max = 0;
    for (i[0] = 0; i[0] < n[0]; i[0]++)
    for (i[1] = 0; i[1] < n[1]; i[1]++)
    for (i[2] = 0; i[2] < n[2]; i[2]++) {
      double r = v[0][i[0]+off[0]] * v[1][i[1]+off[1]] * v[2][i[2]+off[2]];
      double exact = lam == 0 ? 0 : r / lam;
      double e = fabs(phi[idx(i, gp, sp)] - exact);
      /* a NaN must fail the check */
      if (e > err || e != e) err = e;
      if (fabs(r) > rhs_max) rhs_max = fabs(r);
      /* mark the unknown to check the ghost points below */
      phi[idx(i, gp, sp)] = ghost_phi;
    }
    err = global_max(err);
    rhs_max = global_max(rhs_max);

    ghost_ok = 1;
    for (k = 0; k < np; k++)
      if (phi[k] != ghost_phi) ghost_ok = 0;

    sprintf(label, "n=%dx%dx%d mode %d,%d,%d", gn[0], gn[1], gn[2], m[0], m[1], m[2]);
    ut_c_check_err(label, err * lam_min / rhs_max, 100 * REAL_EPS);
    sprintf(label, "n=%dx%dx%d mode %d,%d,%d ghost points of Phi unchanged",
            gn[0], gn[1], gn[2], m[0], m[1], m[2]);
    ut_c_check(global_all(ghost_ok), label);
  }

  solver_finalize(&S);
  ut_c_check(global_all(S.D == NULL), "finalize clears the solver");

  for (d = 0; d < 3; d++) free(v[d]);
  free(phi);
  free(rhs);
}

int main(int argc, char **argv)
{
  const int P = POISFFT_PERIODIC, D = POISFFT_DIRICHLET, N = POISFFT_NEUMANN,
            Ds = POISFFT_DIRICHLET_STAG, Ns = POISFFT_NEUMANN_STAG;
  const int approxs[3] = {POISFFT_SPECTRAL, POISFFT_FINITE_DIFFERENCE_2,
                          POISFFT_FINITE_DIFFERENCE_4};
  /* global sizes and boundary conditions in the C order */
  const int n1[1] = {7};
  const int bcs1[3][2] = {{D, N}, {P, P}, {Ns, Ds}};
  const int n2[2] = {7, 6};
  const int bcs2[3][4] = {{P, P, D, N}, {Ns, Ns, Ds, Ds}, {N, D, P, P}};
  const int n3[3] = {7, 6, 5};
  const int bcs3[4][6] = {{D, D, P, P, N, N}, {Ds, Ns, Ns, Ns, P, P},
                          {P, P, P, P, P, P}, {N, N, Ds, Ds, D, N}};
  int ia, ic, ghosts;
  int rank = 0;
  char name[80];

#ifdef MPI
  int nprocs, npy = 1, npz, provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  /* the solvers call MPI from inside OpenMP single regions */
  if (provided < MPI_THREAD_SERIALIZED) {
    if (rank == 0) fprintf(stderr, "MPI does not provide MPI_THREAD_SERIALIZED\n");
    MPI_Abort(MPI_COMM_WORLD, 5);
  }
  npz = nprocs;
  if (argc >= 3) {
    npy = atoi(argv[1]);
    npz = atoi(argv[2]);
  }
  if (npy * npz != nprocs) {
    if (rank == 0)
      fprintf(stderr, "Process grid %d x %d does not match the number of processes %d\n",
              npy, npz, nprocs);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  /* PFFT expects the process mesh in the C order */
  if (pfft_create_procmesh_2d(MPI_COMM_WORLD, npz, npy, &comm_3d) != 0 ||
      pfft_create_procmesh_1d(MPI_COMM_WORLD, nprocs, &comm_2d) != 0) {
    if (rank == 0) fprintf(stderr, "pfft_create_procmesh failed\n");
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  sprintf(name, "test_c_api " PREC_NAME " grid %dx%d", npy, npz);
#else
  (void) argc;
  (void) argv;
  sprintf(name, "test_c_api " PREC_NAME);
#endif

  ut_c_init(name, rank != 0);

  for (ghosts = 0; ghosts <= 1; ghosts++)
    for (ia = 0; ia < 3; ia++) {
      for (ic = 0; ic < 3; ic++)
        test_problem(1, n1, bcs1[ic], approxs[ia], ghosts);
      for (ic = 0; ic < 3; ic++)
        test_problem(2, n2, bcs2[ic], approxs[ia], ghosts);
      for (ic = 0; ic < 4; ic++)
        test_problem(3, n3, bcs3[ic], approxs[ia], ghosts);
    }

#ifdef MPI
  MPI_Finalize();
#endif
  ut_c_finish();
  return 0;
}
