module PoisFFT
#ifdef MPI
  use PFFT
#endif

  use PoisFFT_Precisions
  
  use PoisFFT_Parameters

  use PoisFFT_SP, PoisFFT_Solver1D_SP => PoisFFT_Solver1D, &
                  PoisFFT_Solver2D_SP => PoisFFT_Solver2D, &
                  PoisFFT_Solver3D_SP => PoisFFT_Solver3D, &
                  PoisFFT_Solver3D_nonuniform_z_SP => PoisFFT_Solver3D_nonuniform_z

  use PoisFFT_DP, PoisFFT_Solver1D_DP => PoisFFT_Solver1D, &
                  PoisFFT_Solver2D_DP => PoisFFT_Solver2D, &
                  PoisFFT_Solver3D_DP => PoisFFT_Solver3D, &
                  PoisFFT_Solver3D_nonuniform_z_DP => PoisFFT_Solver3D_nonuniform_z
contains
#ifdef MPI    
    !Cartesian communicator for the solvers, np(1) processes for the 1D and 2D
    !solvers, np(1:2) for the 3D solvers (in the C order: z, y).
    !ierror /= 0 if the number of processes does not match or on an error.
    subroutine PoisFFT_InitMPIGrid(MPI_comm, np, PoisFFT_comm, ierror)
      use mpi, only: MPI_Comm_size
      integer, intent(in) :: mpi_comm
      integer(c_int), intent(in) :: np(:)
      integer, intent(out) :: PoisFFT_comm, ierror
      integer :: comm_size
      
      call MPI_Comm_size(mpi_comm, comm_size, ierror)
      
      if (ierror/=0 .or. comm_size/=product(np) .or. &
          size(np) < 1 .or. size(np) > 2) then
        ierror = 1
        return
      end if
      
      !This will fail to compile if integer /= 4 bytes as PFFT assumes.
      if (size(np)==1) then
        ierror = pfft_create_procmesh_1d(MPI_comm, &
                                         np(1), &
                                         PoisFFT_comm)
      else
        ierror = pfft_create_procmesh_2d(MPI_comm, &
                                         np(1), &
                                         np(2), &
                                         PoisFFT_comm)
      end if

    end subroutine
    
    !The local block of the PFFT decomposition (PFFT_TRANSPOSED_NONE) of the
    !global grid Nos(1:rnk), rnk = 2 or 3, in the Fortran order. local_ni and
    !local_i_start are the local sizes and offsets for the solvers.
    subroutine PoisFFT_LocalGridSize(rnk,Nos,PoisFFT_comm,local_ni,local_i_start,local_no,local_o_start)
      integer(c_int), value :: rnk
      integer(c_intptr_t), intent(in) :: Nos(:)
      integer, intent(in) :: PoisFFT_comm
      integer(c_intptr_t), intent(out) :: local_ni(:)
      integer(c_intptr_t), intent(out) :: local_i_start(:)
      integer(c_intptr_t), intent(out) :: local_no(:)
      integer(c_intptr_t), intent(out) :: local_o_start(:)
      integer(c_intptr_t) :: alloc_local
      integer(c_intptr_t) :: Nback(1:rnk)
      
      !PFFT uses the C order of the dimensions
      Nback = Nos(rnk:1:-1)
      alloc_local = pfft_local_size_dft(rnk, Nback, PoisFFT_comm, PFFT_TRANSPOSED_NONE, &
                                        local_ni, local_i_start, local_no, local_o_start)
      local_ni(1:rnk) = local_ni(rnk:1:-1)
      local_i_start(1:rnk) = local_i_start(rnk:1:-1)
      local_no(1:rnk) = local_no(rnk:1:-1)
      local_o_start(1:rnk) = local_o_start(rnk:1:-1)

    end subroutine
#endif
end module PoisFFT
