#if dimensions == 1

#define PoisFFT_SolverXD PoisFFT_Solver1D
#define PoisFFT_PlanXD   PoisFFT_Plan1D
#define realplantypes    plantypes(1)
#define nxyzs            D%gnx
#define colons           :

#elif dimensions == 2

#define PoisFFT_SolverXD PoisFFT_Solver2D
#define PoisFFT_PlanXD   PoisFFT_Plan2D
#define realplantypes    plantypes(1),plantypes(2)
#define nxyzs            D%gny,D%gnx
#define colons           :,:

#else

#define PoisFFT_SolverXD PoisFFT_Solver3D
#define PoisFFT_PlanXD   PoisFFT_Plan3D
#define realplantypes    plantypes(1),plantypes(2),plantypes(3)
#define nxyzs            D%gnz,D%gny,D%gnx
#define colons           :,:,:

#endif

#if CP == 4
#define pfft_cmplx pfftf_plan_dft_3d
#define pfft_real pfftf_plan_r2r_3d
#else
#define pfft_cmplx pfft_plan_dft_3d
#define pfft_real pfft_plan_r2r_3d
#endif

      type(PoisFFT_PlanXD) :: plan
      
      type(PoisFFT_SolverXD), intent(inout) :: D
      integer, intent(in), dimension(:)     :: plantypes

      if (plantypes(1)==FFT_Complex) then
       
        if (size(plantypes)<2) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif
       
        plan%dir = plantypes(2)

#if MPI && dimensions >=3
! print *,"forw",int([nxyzs],c_intptr_t), &
!           size(D%cwork), size(D%cwork), D%mpi_comm, &
!           plan%dir, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT
        plan%planptr = pfft_cmplx(int([nxyzs],c_intptr_t), &
          D%cwork, D%cwork, D%mpi_comm, &
          plan%dir, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
#else
        plan%planptr = fftw_plan_gen(nxyzs , D%cwork, D%cwork,&
                        plan%dir, FFTW_MEASURE)
#endif
        
        
      else
      
        if (size(plantypes)< dimensions ) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif


#if MPI && dimensions >=3
        plan%planptr = pfft_real(int([nxyzs],c_intptr_t), &
          D%rwork, D%rwork, D%mpi_comm, &
          plantypes, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
#else
        plan%planptr = fftw_plan_gen(nxyzs , D%rwork, D%rwork,&
                        realplantypes , FFTW_MEASURE)
#endif
        
      endif
      
      plan%planowner=.true.
      

#undef colons
#undef realplantypes
#undef nxyzs
#undef PoisFFT_SolverXD
#undef PoisFFT_PlanXD
#undef SliceXD
#undef pfft_cmplx
#undef pfft_real