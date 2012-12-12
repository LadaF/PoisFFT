#if dimensions == 1

#define PoisFFT_SolverXD PoisFFT_Solver1D
#define PoisFFT_PlanXD   PoisFFT_Plan1D
#define SliceXD          Slice1D
#define realplantypes    plantypes(1)
#define nxyzs            plan%nx
#define colons           :

#elif dimensions == 2

#define PoisFFT_SolverXD PoisFFT_Solver2D
#define PoisFFT_PlanXD   PoisFFT_Plan2D
#define SliceXD          Slice2D
#define realplantypes    plantypes(1),plantypes(2)
#define nxyzs            plan%ny,plan%nx
#define colons           :,:

#else

#define PoisFFT_SolverXD PoisFFT_Solver3D
#define PoisFFT_PlanXD   PoisFFT_Plan3D
#define SliceXD          Slice3D
#define realplantypes    plantypes(1),plantypes(2),plantypes(3)
#define nxyzs            plan%nz,plan%ny,plan%nx
#define colons           :,:,:

#endif



      type(PoisFFT_PlanXD) :: plan
      
      type(PoisFFT_SolverXD), intent(inout) :: D
      type(SliceXD), intent(in)            :: sl
      logical, intent(in)                   :: usework
      integer, intent(in), dimension(:)     :: plantypes

      complex(cprec), dimension(colons),&!contiguous,
                                     pointer :: data_loc_complex
      real(rprec), dimension(colons),&!contiguous,
                                     pointer :: data_loc_real

      nullify(data_loc_complex)
      nullify(data_loc_real)
      
      plan%nx = (sl%ie - sl%is + 1)
      plan%cnt = plan%nx     
#if dimensions >= 2
      plan%ny = (sl%je - sl%js + 1)
      plan%cnt = plan%cnt * plan%ny
#endif
#if dimensions >= 3
      plan%nz = (sl%ke - sl%ks + 1)
      plan%cnt = plan%cnt * plan%nz
#endif

      if (plantypes(1)==FFT_Complex) then
       
        if (size(plantypes)<2) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif
       
        plan%dir = plantypes(2)

        if (usework) then
          data_loc_complex => D%cwork
        else
          call data_allocate_complex(D, data_loc_complex)
        endif


        plan%planptr = fftw_plan_gen(nxyzs , data_loc_complex, data_loc_complex,&
                        plan%dir, FFTW_MEASURE)

        if (.not.usework) then
          call data_deallocate(data_loc_complex)
        endif
        
      else
      
        if (size(plantypes)< dimensions ) then
          write (*,*) "Error: not enough flags when creating PoisFFT_PlanXD"
          STOP
        endif

        if (usework) then
          data_loc_real => D%rwork
        else
          call data_allocate_real(D, data_loc_real)
        endif


        plan%planptr = fftw_plan_gen(nxyzs , data_loc_real, data_loc_real,&
                        realplantypes , FFTW_MEASURE)

        if (.not.usework) then
          call data_deallocate(data_loc_real)
        endif
        
      endif
      
      plan%planowner=.true.
      

#undef colons
#undef realplantypes
#undef nxyzs
#undef PoisFFT_SolverXD
#undef PoisFFT_PlanXD
#undef SliceXD
