!Combinations of boundary conditions of the 3D tests, shared by the serial
!and the MPI tests. Every specialized solver routine of PoisFFT_Solver3D
!(named in the comments) is covered with its variants. Under MPI the
!routines NsPP, PNsP, DDP, 2real1real, 1real2real and Separable use the
!general separable MPI solver. Separable takes every combination without a
!specialized routine, cases 32..52.
module bc_combos
  use PoisFFT_Parameters
  implicit none

  private
  public :: n_combos_3d, combos_3d

  integer, parameter :: P = PoisFFT_Periodic, D = PoisFFT_Dirichlet, &
                        N = PoisFFT_Neumann, Ds = PoisFFT_DirichletStag, &
                        Ns = PoisFFT_NeumannStag

  integer, parameter :: n_combos_3d = 52
  integer, parameter :: combos_3d(6, n_combos_3d) = reshape([ &
      P, P, P, P, P, P,           & ! 1 FullPeriodic
      D, D, D, D, D, D,           & ! 2 FullDirichlet
      Ds,Ds,Ds,Ds,Ds,Ds,          & ! 3
      N, N, N, N, N, N,           & ! 4 FullNeumann
      Ns,Ns,Ns,Ns,Ns,Ns,          & ! 5
      P, P, P, P, N, N,           & ! 6 PPNs
      P, P, P, P, Ns,Ns,          & ! 7
      P, P, P, P, D, D,           & ! 8 PPD
      P, P, P, P, Ds,Ds,          & ! 9
      P, P, P, P, Ds,Ns,          & !10
      P, P, P, P, Ns,Ds,          & !11
      P, P, P, P, D, N,           & !12
      P, P, P, P, N, D,           & !13
      N, N, P, P, P, P,           & !14 NsPP
      Ns,Ns,P, P, P, P,           & !15
      P, P, N, N, P, P,           & !16 PNsP
      P, P, Ns,Ns,P, P,           & !17
      P, P, N, N, N, N,           & !18 PNsNs
      P, P, Ns,Ns,Ns,Ns,          & !19
      D, D, D, D, P, P,           & !20 DDP
      Ds,Ds,Ds,Ds,P, P,           & !21
      D, D, D, D, N, N,           & !22 2real1real
      Ds,Ds,Ds,Ds,Ns,Ns,          & !23
      Ns,Ns,Ns,Ns,Ds,Ds,          & !24
      Ds,Ns,Ds,Ns,Ns,Ns,          & !25
      Ds,Ds,Ds,Ds,Ns,Ds,          & !26
      N, N, D, D, D, D,           & !27 1real2real
      Ns,Ns,Ds,Ds,Ds,Ds,          & !28
      Ds,Ds,Ns,Ns,Ns,Ns,          & !29
      Ns,Ds,Ds,Ds,Ds,Ds,          & !30
      D, D, N, N, N, N,           & !31
      D, D, P, P, D, D,           & !32 Separable: real, periodic, real
      Ds,Ds,P, P, Ds,Ds,          & !33
      N, N, P, P, N, N,           & !34
      Ns,Ns,P, P, Ns,Ns,          & !35
      D, D, P, P, N, N,           & !36
      Ns,Ns,P, P, Ds,Ds,          & !37
      N, D, P, P, D, N,           & !38
      Ds,Ns,P, P, Ns,Ds,          & !39
      D, D, P, P, P, P,           & !40 Separable: real, periodic, periodic
      N, D, P, P, P, P,           & !41
      P, P, Ds,Ds,P, P,           & !42 Separable: periodic, real, periodic
      P, P, D, N, P, P,           & !43
      P, P, D, D, N, N,           & !44 Separable: periodic, real, real
      P, P, Ns,Ns,Ds,Ns,          & !45
      N, N, N, N, P, P,           & !46 Separable: real, real, periodic
      D, D, Ns,Ds,P, P,           & !47
      D, D, N, N, D, D,           & !48 Separable: real in all directions
      D, D, Ns,Ns,N, D,           & !49
      N, N, N, N, Ns,Ns,          & !50
      N, D, N, D, N, D,           & !51
      Ns,Ns,N, N, Ns,Ns           & !52
    ], [6, n_combos_3d])

end module bc_combos
