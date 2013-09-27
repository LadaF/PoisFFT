module FFT_SP
#define RP 4
#define CP 4

#include "fft-inc.f90"

#undef RP
#undef CP
end module FFT_SP

module FFT_DP
#define RP 8
#define CP 8

#include "fft-inc.f90"

#undef RP
#undef CP
end module FFT_DP
