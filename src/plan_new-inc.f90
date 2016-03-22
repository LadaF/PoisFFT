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

#if (PREC == 2)
#define pfft_cmplx pfft_plan_dft
#define pfft_real pfft_plan_r2r
#define fftw_cmplx_mpi_2d fftw_mpi_plan_dft_2d
#else
#define pfft_cmplx pfftf_plan_dft
#define pfft_real pfftf_plan_r2r
#define fftw_cmplx_mpi_2d fftwf_mpi_plan_dft_2d
#endif

      type(PoisFFT_PlanXD) :: plan
      
      type(PoisFFT_SolverXD), intent(inout) :: D
      integer, intent(in), dimension(:)     :: plantypes
      logical, intent(in), optional         :: distributed
      logical :: distr

      distr = .false.

      if (plantypes(1)==FFT_Complex) then

        if (size(plantypes)<2) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif
       
        plan%dir = plantypes(2)

#if defined(MPI) && dimensions > 1
        if (present(distributed)) then
          distr = distributed
        else
          distr = .true.
        end if

#if dimensions == 2
        if (distr) then
          if (D%mpi%comm_dim==1) then
            plan%planptr = fftw_cmplx_mpi_2d(int(D%gny,c_intptr_t),int(D%gnx,c_intptr_t), &
              D%cwork, D%cwork, D%mpi%comm, &
              plan%dir, FFTW_MEASURE)
            plan%method = FFT_DISTRIBUTED_FFTW
          else
            plan%planptr = pfft_cmplx(dimensions,int([nxyzs],c_intptr_t), &
              D%cwork, D%cwork, D%mpi%comm, &
              plan%dir, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
          end if
#elif dimensions == 3
        if (distr) then
          plan%planptr = pfft_cmplx(dimensions,int([nxyzs],c_intptr_t), &
            D%cwork, D%cwork, D%mpi%comm, &
            plan%dir, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
#endif
        else
          plan%planptr = fftw_plan_gen(nxyzs , D%cwork, D%cwork,&
                          plan%dir, FFTW_MEASURE)
        end if
#else
        plan%planptr = fftw_plan_gen(nxyzs , D%cwork, D%cwork,&
                        plan%dir, FFTW_MEASURE)
#endif

      else

        if (size(plantypes)< dimensions ) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD, there must be one per dimension."
          STOP
        endif


#if defined(MPI) && dimensions > 1
        if (present(distributed)) then
          distr = distributed
        else
          distr = .true.
        end if

        if (distr) then
          plan%planptr = pfft_real(dimensions,int([nxyzs],c_intptr_t), &
            D%rwork, D%rwork, D%mpi%comm, &
            plantypes, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
        else
          plan%planptr = fftw_plan_gen(nxyzs , D%rwork, D%rwork,&
                          realplantypes , FFTW_MEASURE)
        end if
#else
        plan%planptr = fftw_plan_gen(nxyzs , D%rwork, D%rwork,&
                        realplantypes , FFTW_MEASURE)
#endif

      endif
      
      plan%planowner = .true.

      plan%distributed = distr

      if (.not.c_associated(plan%planptr)) stop "Error, FFT plan not created!"

#undef colons
#undef realplantypes
#undef nxyzs
#undef PoisFFT_SolverXD
#undef PoisFFT_PlanXD
#undef SliceXD
#undef pfft_cmplx
#undef pfft_real
#undef fftw_cmplx_mpi_2d