module FFT_SP
#define RP SRP
#define CP SCP

#include "fft-inc.f90"

#undef RP
#undef CP
end module FFT_SP

module FFT_DP
#define RP DRP
#define CP DCP

#include "fft-inc.f90"

#undef RP
#undef CP
end module FFT_DP
