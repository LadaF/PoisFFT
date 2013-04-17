module PoisFFT

  use PoisFFT_Precisions
  
  use PoisFFT_Parameters

  use PoisFFT_SP, PoisFFT_Solver1D_SP => PoisFFT_Solver1D, &
                  PoisFFT_Solver2D_SP => PoisFFT_Solver2D, &
                  PoisFFT_Solver3D_SP => PoisFFT_Solver3D

  use PoisFFT_DP, PoisFFT_Solver1D_DP => PoisFFT_Solver1D, &
                  PoisFFT_Solver2D_DP => PoisFFT_Solver2D, &
                  PoisFFT_Solver3D_DP => PoisFFT_Solver3D


end module PoisFFT