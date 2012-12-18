
module PoisFFT_SP
  use FFT_SP
#define RP SRP
#define CP SCP

#define _RP _SRP
#define _CP _SCP

#include "poisfft-inc.f90"

#undef RP
#undef CP
#undef _RP
#undef _CP
end module PoisFFT_SP

module PoisFFT_DP
  use FFT_DP
#define RP DRP
#define CP DCP
#define _RP _DRP
#define _CP _DCP

#include "poisfft-inc.f90"

#undef RP
#undef CP
#undef _RP
#undef _CP
end module PoisFFT_DP
