module FFT_SP
#define PREC 1

#include "fft-inc.f90"

#undef PREC
end module FFT_SP

module FFT_DP
#define PREC 2

#include "fft-inc.f90"

#undef PREC
end module FFT_DP
