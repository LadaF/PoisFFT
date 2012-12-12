  integer, parameter :: C_FFTW_R2R_KIND = C_INT32_T

  integer(C_INT), parameter :: FFTW_R2HC = 0
  integer(C_INT), parameter :: FFTW_HC2R = 1
  integer(C_INT), parameter :: FFTW_DHT = 2
  integer(C_INT), parameter :: FFTW_REDFT00 = 3
  integer(C_INT), parameter :: FFTW_REDFT01 = 4
  integer(C_INT), parameter :: FFTW_REDFT10 = 5
  integer(C_INT), parameter :: FFTW_REDFT11 = 6
  integer(C_INT), parameter :: FFTW_RODFT00 = 7
  integer(C_INT), parameter :: FFTW_RODFT01 = 8
  integer(C_INT), parameter :: FFTW_RODFT10 = 9
  integer(C_INT), parameter :: FFTW_RODFT11 = 10
  integer(C_INT), parameter :: FFTW_FORWARD = -1
  integer(C_INT), parameter :: FFTW_BACKWARD = +1
  integer(C_INT), parameter :: FFTW_MEASURE = 0
  integer(C_INT), parameter :: FFTW_DESTROY_INPUT = 1
  integer(C_INT), parameter :: FFTW_UNALIGNED = 2
  integer(C_INT), parameter :: FFTW_CONSERVE_MEMORY = 4
  integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8
  integer(C_INT), parameter :: FFTW_PRESERVE_INPUT = 16
  integer(C_INT), parameter :: FFTW_PATIENT = 32
  integer(C_INT), parameter :: FFTW_ESTIMATE = 64
  integer(C_INT), parameter :: FFTW_ESTIMATE_PATIENT = 128
  integer(C_INT), parameter :: FFTW_BELIEVE_PCOST = 256
  integer(C_INT), parameter :: FFTW_NO_DFT_R2HC = 512
  integer(C_INT), parameter :: FFTW_NO_NONTHREADED = 1024
  integer(C_INT), parameter :: FFTW_NO_BUFFERING = 2048
  integer(C_INT), parameter :: FFTW_NO_INDIRECT_OP = 4096
  integer(C_INT), parameter :: FFTW_ALLOW_LARGE_GENERIC = 8192
  integer(C_INT), parameter :: FFTW_NO_RANK_SPLITS = 16384
  integer(C_INT), parameter :: FFTW_NO_VRANK_SPLITS = 32768
  integer(C_INT), parameter :: FFTW_NO_VRECURSE = 65536
  integer(C_INT), parameter :: FFTW_NO_SIMD = 131072
  integer(C_INT), parameter :: FFTW_NO_SLOW = 262144
  integer(C_INT), parameter :: FFTW_NO_FIXED_RADIX_LARGE_N = 524288
  integer(C_INT), parameter :: FFTW_ALLOW_PRUNING = 1048576
  integer(C_INT), parameter :: FFTW_WISDOM_ONLY = 2097152


  
  interface fftw_execute_gen
    subroutine fftw_execute_r2r(p,in,out) bind(C, name='fftw_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine fftw_execute_r2r
    subroutine fftwf_execute_r2r(p,in,out) bind(C, name='fftwf_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine fftwf_execute_r2r
    subroutine fftw_execute_dft(p,in,out) bind(C, name='fftw_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftw_execute_dft
    subroutine fftwf_execute_dft(p,in,out) bind(C, name='fftwf_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftwf_execute_dft
  end interface fftw_execute_gen


  

  interface fftw_plan_gen
    type(C_PTR) function fftw_plan_dft_1d(n,in,out,sign,flags) bind(C, name='fftw_plan_dft_1d')
      import
      integer(C_INT), value :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_dft_1d

    type(C_PTR) function fftw_plan_dft_2d(n0,n1,in,out,sign,flags) bind(C, name='fftw_plan_dft_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(n1,*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(n1,*), intent(out) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_dft_2d

    type(C_PTR) function fftw_plan_dft_3d(n0,n1,n2,in,out,sign,flags) bind(C, name='fftw_plan_dft_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(n2,n1,*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(n2,n1,*), intent(out) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_dft_3d

    type(C_PTR) function fftwf_plan_dft_1d(n,in,out,sign,flags) bind(C, name='fftwf_plan_dft_1d')
      import
      integer(C_INT), value :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_1d

    type(C_PTR) function fftwf_plan_dft_2d(n0,n1,in,out,sign,flags) bind(C, name='fftwf_plan_dft_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(n1,*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(n1,*), intent(out) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_2d

    type(C_PTR) function fftwf_plan_dft_3d(n0,n1,n2,in,out,sign,flags) bind(C, name='fftwf_plan_dft_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(n2,n1,*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(n2,n1,*), intent(out) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_3d

    type(C_PTR) function fftw_plan_r2r_1d(n,in,out,kind,flags) bind(C, name='fftw_plan_r2r_1d')
      import
      integer(C_INT), value :: n
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_FFTW_R2R_KIND), value :: kind
      integer(C_INT), value :: flags
    end function fftw_plan_r2r_1d

    type(C_PTR) function fftw_plan_r2r_2d(n0,n1,in,out,kind0,kind1,flags) bind(C, name='fftw_plan_r2r_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      real(C_DOUBLE), dimension(n1,*), intent(out) :: in
      real(C_DOUBLE), dimension(n1,*), intent(out) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftw_plan_r2r_2d

    type(C_PTR) function fftw_plan_r2r_3d(n0,n1,n2,in,out,kind0,kind1,kind2,flags) bind(C, name='fftw_plan_r2r_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      real(C_DOUBLE), dimension(n2,n1,*), intent(out) :: in
      real(C_DOUBLE), dimension(n2,n1,*), intent(out) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftw_plan_r2r_3d

    type(C_PTR) function fftwf_plan_r2r_1d(n,in,out,kind,flags) bind(C, name='fftwf_plan_r2r_1d')
      import
      integer(C_INT), value :: n
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_FFTW_R2R_KIND), value :: kind
      integer(C_INT), value :: flags
    end function fftwf_plan_r2r_1d

    type(C_PTR) function fftwf_plan_r2r_2d(n0,n1,in,out,kind0,kind1,flags) bind(C, name='fftwf_plan_r2r_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      real(C_FLOAT), dimension(n1,*), intent(out) :: in
      real(C_FLOAT), dimension(n1,*), intent(out) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftwf_plan_r2r_2d

    type(C_PTR) function fftwf_plan_r2r_3d(n0,n1,n2,in,out,kind0,kind1,kind2,flags) bind(C, name='fftwf_plan_r2r_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      real(C_FLOAT), dimension(n2,n1,*), intent(out) :: in
      real(C_FLOAT), dimension(n2,n1,*), intent(out) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftwf_plan_r2r_3d

  end interface fftw_plan_gen


  interface

    type(C_PTR) function fftw_malloc(n) bind(C, name='fftw_malloc')
      import
      integer(C_SIZE_T), value :: n
    end function fftw_malloc

    subroutine fftw_free(p) bind(C, name='fftw_free')
      import
      type(C_PTR), value :: p
    end subroutine fftw_free

    subroutine fftwf_destroy_plan(p) bind(C, name='fftwf_destroy_plan')
      import
      type(C_PTR), value :: p
    end subroutine fftwf_destroy_plan

    subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
      import
      type(C_PTR), value :: p
    end subroutine fftw_destroy_plan


    subroutine fftw_plan_with_nthreads(nthreads) bind(C, name='fftw_plan_with_nthreads')
      import
      integer(C_INT), value :: nthreads
    end subroutine fftw_plan_with_nthreads

    integer(C_INT) function fftw_init_threads() bind(C, name='fftw_init_threads')
      import
    end function fftw_init_threads

    subroutine fftw_cleanup_threads() bind(C, name='fftw_cleanup_threads')
      import
    end subroutine fftw_cleanup_threads

   end interface
  