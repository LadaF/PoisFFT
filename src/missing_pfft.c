#if defined(MPI) && defined(MISSING_PFFT_R2R)
#include "pfft.h"

#define PFFT_DEFINE_MISSING(PX, X, R, C, INT)                                               \
PX(plan) PX(plan_r2r_f03)(int rnk, const INT * Nos, R * in, R * out, MPI_Fint f_comm_cart, const PX(r2r_kind) * kinds, unsigned pfft_flags) \
{\
  MPI_Comm comm_cart;\
\
  comm_cart = MPI_Comm_f2c(f_comm_cart);\
  PX(plan) ret = PX(plan_r2r)(rnk, Nos, in, out, comm_cart, kinds, pfft_flags);\
  return ret;\
}

PFFT_DEFINE_MISSING(PFFT_MANGLE_DOUBLE, FFTW_MANGLE_DOUBLE, double, pfft_complex, ptrdiff_t)
PFFT_DEFINE_MISSING(PFFT_MANGLE_FLOAT, FFTW_MANGLE_FLOAT, float, pfftf_complex, ptrdiff_t)
#endif