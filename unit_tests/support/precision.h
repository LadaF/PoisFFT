!Selects the precision of a test program. Compile with -DSINGLE for the
!single precision variant, the default is double precision.
#ifdef SINGLE
#define SOLVER1D PoisFFT_Solver1D_SP
#define SOLVER2D PoisFFT_Solver2D_SP
#define SOLVER3D PoisFFT_Solver3D_SP
#define SOLVER3D_NUZ PoisFFT_Solver3D_nonuniform_z_SP
#define FFT_MODULE FFT_SP
#define RP_KIND SRP
#define CP_KIND SCP
#define PREC_NAME "SP"
#else
#define SOLVER1D PoisFFT_Solver1D_DP
#define SOLVER2D PoisFFT_Solver2D_DP
#define SOLVER3D PoisFFT_Solver3D_DP
#define SOLVER3D_NUZ PoisFFT_Solver3D_nonuniform_z_DP
#define FFT_MODULE FFT_DP
#define RP_KIND DRP
#define CP_KIND DCP
#define PREC_NAME "DP"
#endif
