module PoisFFT_Parameters

  integer, parameter :: PoisFFT_Periodic = 0
  integer, parameter :: PoisFFT_Dirichlet = 1
  integer, parameter :: PoisFFT_Neumann = 2
  integer, parameter :: PoisFFT_DirichletStag = 3
  integer, parameter :: PoisFFT_NeumannStag = 4

  integer, parameter :: PoisFFT_Spectral = 0
  integer, parameter :: PoisFFT_FiniteDifference2 = 2
  integer, parameter :: PoisFFT_FiniteDifference4 = 4

  !planning of the FFTW and PFFT plans: the measured fastest plans
  !(FFTW_MEASURE) can differ from run to run and with them the rounding of
  !the results, the estimated plans (FFTW_ESTIMATE) are reproducible
  integer, parameter :: PoisFFT_Measure = 0
  integer, parameter :: PoisFFT_Estimate = 1
end module PoisFFT_Parameters