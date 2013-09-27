#if dimensions == 1
#define PoisFFT_SolverXD PoisFFT_Solver1D
#define dims [D%nx]
#define colons :
#elif dimensions == 2
#define PoisFFT_SolverXD PoisFFT_Solver2D
#define dims [D%nx,D%ny]
#define colons :,:
#else
#define PoisFFT_SolverXD PoisFFT_Solver3D
#define dims [D%nx,D%ny,D%nz]
#define colons :,:,:
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
!       datatype , dimension( colons ), pointer, optional    :: data
      type(c_ptr) :: p
#if MPI && dimensions == 3
      integer(c_size_t) :: cnt
      integer(c_intptr_t) :: a1(3),a2(3),a3(3),a4(3)
      cnt =  pfft_local_size_dft_3d(int([D%gnz,D%gny,D%gnx],c_intptr_t), &
                  D%mpi_comm, PFFT_TRANSPOSED_NONE, &
                  a1,a2,a3,a4)
      if (.not.all(a1(dimensions:1:-1)==dims)) stop "Error. Inconsistent size of local arrays!"
      p = fftw_malloc( sizeof( dataf ) * cnt )
#else
      p = fftw_malloc( sizeof( dataf ) * D%cnt )
#endif

      if (c_associated(p)) then

          call c_f_pointer(p, D%xwork, dims )

      else
        stop "Data allocate error, fftw_malloc returned NULL."
      endif

      D%xwork = 0
#undef xwork
#undef datatype
#undef dataf
#undef colons
#undef dims
#undef PoisFFT_SolverXD
