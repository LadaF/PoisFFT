
module PoisFFT_SP
  use FFT_SP
#define PREC 1

#include "poisfft-inc.f90"

#undef PREC
end module PoisFFT_SP

module PoisFFT_DP
  use FFT_DP
#define PREC 2

#include "poisfft-inc.f90"

#undef PREC
end module PoisFFT_DP
