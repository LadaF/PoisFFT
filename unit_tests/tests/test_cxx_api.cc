/* Unit test of the C++ interface, the class template PoisFFT::Solver of
   poisfft.h. The C functions under it are tested by test_c_api.c, this
   test checks the class for 1D, 2D and 3D: the constructor with the
   default arguments and with an approximation, execute() with and
   without ghost points, the destructor at the end of the scope and that
   the class cannot be copied.
   Compiled with -DSINGLE for float. See test_c_api.c for the method. */
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <type_traits>
#include <vector>
#include "poisfft.h"
#include "c_test_support.h"

#ifdef SINGLE
typedef float real;
#define PREC_NAME "SP"
#else
typedef double real;
#define PREC_NAME "DP"
#endif

namespace {

// the solvers own the Fortran object and must not be copied
static_assert(!std::is_copy_constructible<PoisFFT::Solver<3, real> >::value &&
              !std::is_copy_assignable<PoisFFT::Solver<3, real> >::value,
              "PoisFFT::Solver must not be copyable");

const double Ls[3] = {2.3, 1.7, 1.1};
const int ng_phi[3] = {2, 1, 0}, ng_rhs[3] = {1, 0, 2};
const real ghost_rhs = -777, ghost_phi = 12345;

/* All modes of one problem, see test_problem in test_c_api.c. */
template <unsigned int dims>
void test_problem(const int *ns, const int *bcs, int approx, bool ghosts,
                  const char *section)
{
  int n[3] = {1, 1, 1}, gp[3] = {0, 0, 0}, gr[3] = {0, 0, 0}, sp[3], sr[3];
  real Lr[dims];
  for (unsigned d = 0; d < dims; d++) {
    int pd = 3 - dims + d;
    n[pd] = ns[d];
    if (ghosts) {
      gp[pd] = ng_phi[d];
      gr[pd] = ng_rhs[d];
    }
    Lr[d] = real(Ls[d]);
  }
  for (int d = 0; d < 3; d++) {
    sp[d] = n[d] + 2 * gp[d];
    sr[d] = n[d] + 2 * gr[d];
  }
  auto idx = [](const int *i, const int *ng, const int *s) {
    return (size_t(i[0] + ng[0]) * s[1] + (i[1] + ng[1])) * s[2] + (i[2] + ng[2]);
  };

  std::vector<real> phi(size_t(sp[0]) * sp[1] * sp[2]), rhs(size_t(sr[0]) * sr[1] * sr[2]);
  std::vector<double> v[3];
  double l[3] = {0, 0, 0};
  for (int d = 0; d < 3; d++) v[d].assign(n[d], 1.0);

  ut_c_section(section);
  double lam_min = ut_c_min_eigenvalue(dims, bcs, ns, Ls, approx);

  auto run = [&](PoisFFT::Solver<dims, real> &S) {
    int m[3], i[3];
    for (m[0] = 1; m[0] <= n[0]; m[0]++)
    for (m[1] = 1; m[1] <= n[1]; m[1]++)
    for (m[2] = 1; m[2] <= n[2]; m[2]++) {
      for (unsigned d = 0; d < dims; d++) {
        int pd = 3 - dims + d;
        ut_c_mode_1d(&bcs[2*d], ns[d], Ls[d], m[pd], approx, v[pd].data(), &l[pd]);
      }
      double lam = l[0] + l[1] + l[2];

      std::fill(rhs.begin(), rhs.end(), ghost_rhs);
      std::fill(phi.begin(), phi.end(), ghost_phi);
      for (i[0] = 0; i[0] < n[0]; i[0]++)
      for (i[1] = 0; i[1] < n[1]; i[1]++)
      for (i[2] = 0; i[2] < n[2]; i[2]++)
        rhs[idx(i, gr, sr)] = real(v[0][i[0]] * v[1][i[1]] * v[2][i[2]]);

      if (ghosts)
        S.execute(phi.data(), rhs.data(), ng_phi, ng_rhs);
      else
        S.execute(phi.data(), rhs.data());

      double err = 0, rhs_max = 0;
      for (i[0] = 0; i[0] < n[0]; i[0]++)
      for (i[1] = 0; i[1] < n[1]; i[1]++)
      for (i[2] = 0; i[2] < n[2]; i[2]++) {
        double r = v[0][i[0]] * v[1][i[1]] * v[2][i[2]];
        double exact = lam == 0 ? 0 : r / lam;
        double e = std::fabs(phi[idx(i, gp, sp)] - exact);
        if (e > err || std::isnan(e)) err = e;
        rhs_max = std::max(rhs_max, std::fabs(r));
        phi[idx(i, gp, sp)] = ghost_phi;
      }
      bool ghost_ok = true;
      for (real x : phi)
        if (x != ghost_phi) ghost_ok = false;

      char label[80];
      std::snprintf(label, sizeof label, "n=%dx%dx%d mode %d,%d,%d", n[0], n[1], n[2], m[0], m[1], m[2]);
      ut_c_check_err(label, err * lam_min / rhs_max, 100 * std::numeric_limits<real>::epsilon());
      std::snprintf(label, sizeof label, "n=%dx%dx%d mode %d,%d,%d ghost points of Phi unchanged",
                    n[0], n[1], n[2], m[0], m[1], m[2]);
      ut_c_check(ghost_ok, label);
    }
  };

  // the destructors at the end of the blocks finalize the solvers
  if (approx == PoisFFT::SPECTRAL) {
    // the default arguments
    PoisFFT::Solver<dims, real> S(ns, Lr, bcs);
    run(S);
  } else {
    PoisFFT::Solver<dims, real> S(ns, Lr, bcs, approx);
    run(S);
  }
}

} // namespace

int main()
{
  using namespace PoisFFT;
  const int n1[1] = {6};
  const int bcs1[2] = {NEUMANN_STAG, DIRICHLET_STAG};
  const int n2[2] = {4, 7};
  const int bcs2[4] = {DIRICHLET, NEUMANN, PERIODIC, PERIODIC};
  const int n3[3] = {5, 3, 4};
  const int bcs3[6] = {PERIODIC, PERIODIC, NEUMANN, NEUMANN, DIRICHLET_STAG, DIRICHLET_STAG};

  ut_c_init("test_cxx_api " PREC_NAME, 0);

  for (int ghosts = 0; ghosts <= 1; ghosts++) {
    test_problem<1>(n1, bcs1, SPECTRAL, ghosts, ghosts ? "C++ 1D Ns,Ds spectral ghosts" : "C++ 1D Ns,Ds spectral");
    test_problem<2>(n2, bcs2, FINITE_DIFFERENCE_2, ghosts, ghosts ? "C++ 2D D,N P,P FD2 ghosts" : "C++ 2D D,N P,P FD2");
    test_problem<3>(n3, bcs3, FINITE_DIFFERENCE_4, ghosts, ghosts ? "C++ 3D P,P N,N Ds,Ds FD4 ghosts" : "C++ 3D P,P N,N Ds,Ds FD4");
    test_problem<3>(n3, bcs3, SPECTRAL, ghosts, ghosts ? "C++ 3D P,P N,N Ds,Ds spectral ghosts" : "C++ 3D P,P N,N Ds,Ds spectral");
  }

  ut_c_finish();
  return 0;
}
