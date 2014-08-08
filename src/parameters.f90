module PoisFFT_Parameters

  integer, parameter :: PoisFFT_Periodic = 0
  integer, parameter :: PoisFFT_Dirichlet = 1
  integer, parameter :: PoisFFT_Neumann = 2
  integer, parameter :: PoisFFT_DirichletStag = 3
  integer, parameter :: PoisFFT_NeumannStag = 4

  integer, parameter :: PoisFFT_Spectral = 0
  integer, parameter :: PoisFFT_FiniteDifference2 = 2
end module PoisFFT_Parameters