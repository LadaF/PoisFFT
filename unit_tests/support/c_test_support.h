/* C and C++ declarations of the test support in c_test_support.f90. */
#ifndef C_TEST_SUPPORT_H
#define C_TEST_SUPPORT_H

#ifdef __cplusplus
extern "C" {
#endif

/* eigenvector v[0..n-1] and eigenvalue of mode m (1..n), 0 for an unsupported pair */
int ut_c_mode_1d(const int *bcs, int n, double L, int m, int approx,
                 double *v, double *lam);
/* smallest non-zero |eigenvalue|, bcs[2*dims], nxyz[dims], Ls[dims] */
double ut_c_min_eigenvalue(int dims, const int *bcs, const int *nxyz,
                           const double *Ls, int approx);

/* silent != 0 suppresses all output (the non-root MPI processes) */
void ut_c_init(const char *name, int silent);
void ut_c_section(const char *name);
void ut_c_check(int condition, const char *label);
void ut_c_check_err(const char *label, double err, double tol);
/* prints the summary, exits with a non-zero code if any check failed */
void ut_c_finish(void);

#ifdef __cplusplus
}
#endif

#endif
