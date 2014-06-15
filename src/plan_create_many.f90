#if dimensions == 1

#define PoisFFT_SolverXD PoisFFT_Solver1D_Many
#define PoisFFT_PlanXD   PoisFFT_Plan1D_Many
#define nxyzs            [D%gnx]
#define colons           :

#elif dimensions == 2

#define PoisFFT_SolverXD PoisFFT_Solver2D_Many
#define PoisFFT_PlanXD   PoisFFT_Plan2D_Many
#define nxyzs            [D%gny,D%gnx]
#define colons           :,:

#endif

#if (PREC == 2)
#define fftw_cmplx fftw_plan_many_dft
#define fftw_real fftw_plan_many_r2r
#define pfft_cmplx pfft_plan_many_dft
#define pfft_real pfft_plan_many_r2r
#else
#define fftw_cmplx fftwf_plan_many_dft
#define fftw_real fftwf_plan_many_r2r
#define pfft_cmplx pfftf_plan_many_dft
#define pfft_real pfftf_plan_many_r2r
#endif

      type(PoisFFT_PlanXD) :: plan
      
      type(PoisFFT_SolverXD), intent(inout) :: D
      integer, intent(in), dimension(:)     :: plantypes
#if defined(MPI) && dimensions > 1
      integer :: i
      integer(c_intptr_t), parameter :: default_blocks(dimensions) = [(0,i = 1, dimensions)]
#endif

      if (plantypes(1)==FFT_Complex) then
       
        if (size(plantypes)<2) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif
       
        plan%dir = plantypes(2)

#if defined(MPI) && dimensions > 1
        plan%planptr = pfft_cmplx(dimensions,int(nxyzs,c_intptr_t),&
                         int(nxyzs,c_intptr_t),int(nxyzs,c_intptr_t), &
                         int(D%howmany,c_intptr_t), default_blocks, default_blocks, &
                         D%cwork, D%cwork, D%mpi%comm, &
                         plan%dir, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
#else
        plan%planptr = fftw_cmplx(dimensions, nxyzs , D%howmany, &
                         D%cwork, nxyzs, 1_c_int, product(nxyzs), &
                         D%cwork, nxyzs, 1_c_int, product(nxyzs), &
                         plan%dir, FFTW_MEASURE)
#endif
        
        
      else
      
        if (size(plantypes)< dimensions ) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif


#if defined(MPI) && dimensions > 1
        plan%planptr = pfft_real(dimensions,int(nxyzs,c_intptr_t),&
                         int(nxyzs,c_intptr_t),int(nxyzs,c_intptr_t), &
                         int(D%howmany,c_intptr_t), default_blocks, default_blocks, &
                         D%rwork, D%rwork, D%mpi%comm, &
                         plantypes, PFFT_TRANSPOSED_NONE + PFFT_MEASURE + PFFT_PRESERVE_INPUT)
#elif defined(MPI) && dimensions == 1
        plan%planptr = fftw_real(dimensions, nxyzs , D%howmany, &
                         D%rwork, nxyzs, 1_c_int, product(nxyzs), &
                         D%rwork, nxyzs, 1_c_int, product(nxyzs), &
                         plantypes, FFTW_MEASURE)
#else
        plan%planptr = fftw_real(dimensions, nxyzs , D%howmany, &
                         D%rwork, nxyzs, 1_c_int, product(nxyzs), &
                         D%rwork, nxyzs, 1_c_int, product(nxyzs), &
                         plantypes, FFTW_MEASURE)
#endif
        
      endif
      
      plan%planowner=.true.
      
      if (.not.c_associated(plan%planptr)) stop "Error, FFT plan not created!"

#undef colons
#undef nxyzs
#undef PoisFFT_SolverXD
#undef PoisFFT_PlanXD
#undef SliceXD
#undef pfft_cmplx
#undef pfft_real
#undef fftw_cmplx
#undef fftw_real
