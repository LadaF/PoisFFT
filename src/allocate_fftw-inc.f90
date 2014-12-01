#if dimensions == 1
#define PoisFFT_SolverXD PoisFFT_Solver1D
#define gdims [D%gnx]
#elif dimensions == 2
#define PoisFFT_SolverXD PoisFFT_Solver2D
#define gdims [D%gny,D%gnx]
#else
#define PoisFFT_SolverXD PoisFFT_Solver3D
#define gdims [D%gnz,D%gny,D%gnx]
#endif

#if realcomplex == 1
#define datatype real(RP)
#define dataf    real(1,RP)
#define xwork    rwork
#else
#define datatype complex(CP)
#define dataf    cmplx(1,1,CP)
#define xwork    cwork
#endif


      type( PoisFFT_SolverXD ), intent(inout)        :: D
      type(c_ptr) :: p
#if defined(MPI) && dimensions>1
      integer(c_size_t) :: cnt
      integer(c_intptr_t) :: a1(dimensions),a2(dimensions),a3(dimensions),a4(dimensions)

      cnt =  pfft_local_size_dft(dimensions,int(gdims,c_intptr_t), &
                  D%mpi%comm, PFFT_TRANSPOSED_NONE, &
                  a1,a2,a3,a4)
      if (.not.all(a1(dimensions:1:-1)==D%nxyz)) stop "Error. Inconsistent size of local arrays!"
      p = fftw_malloc( storage_size( dataf )/storage_size('a') * cnt )
#elif defined(MPI) && dimensions==1
      p = fftw_malloc( storage_size( dataf )/storage_size('a') * int(D%gnx, c_size_t ))
#else
      p = fftw_malloc( storage_size( dataf )/storage_size('a') * D%cnt )
#endif

      if (c_associated(p)) then
#if defined(MPI) && dimensions==1
          call c_f_pointer(p, D%xwork, gdims )
#else
          call c_f_pointer(p, D%xwork, D%nxyz )
#endif
      else
        stop "Data allocate error, fftw_malloc returned NULL."
      endif

      D%xwork = 0
#undef xwork
#undef datatype
#undef dataf
#undef gdims
#undef PoisFFT_SolverXD
