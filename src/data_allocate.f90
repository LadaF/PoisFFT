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
#define datatype real(rprec)
#define dataf    real(1,rprec)
#define xwork    rwork
#else
#define datatype complex(cprec)
#define dataf    cmplx(1,1,cprec)
#define xwork    cwork
#endif


      type( PoisFFT_SolverXD ), intent(inout)        :: D
      datatype , dimension( colons ), pointer, optional    :: data
      type(c_ptr) :: p

      p = fftw_malloc( int(sizeof( dataf ), c_size_t) * int(D%cnt, c_size_t) )

      if (c_associated(p)) then
        if (present(data)) then
          call c_f_pointer(p, data, dims )
        else
          call c_f_pointer(p, D%xwork, dims )
        endif
      else
        stop
      endif
      