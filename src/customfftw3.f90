module fftw3
use iso_c_binding

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
      real(C_DOUBLE), dimension(*), intent(inout) :: out
    end subroutine fftw_execute_r2r
    subroutine fftwf_execute_r2r(p,in,out) bind(C, name='fftwf_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
    end subroutine fftwf_execute_r2r
    subroutine fftw_execute_dft(p,in,out) bind(C, name='fftw_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
    end subroutine fftw_execute_dft
    subroutine fftwf_execute_dft(p,in,out) bind(C, name='fftwf_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
    end subroutine fftwf_execute_dft
    subroutine fftw_execute_dft_r2c(p,in,out) bind(C, name='fftw_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftw_execute_dft_r2c
    subroutine fftw_execute_dft_c2r(p,in,out) bind(C, name='fftw_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine fftw_execute_dft_c2r
    subroutine fftwf_execute_dft_r2c(p,in,out) bind(C, name='fftwf_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
    end subroutine fftwf_execute_dft_r2c
    subroutine fftwf_execute_dft_c2r(p,in,out) bind(C, name='fftwf_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine fftwf_execute_dft_c2r
  end interface fftw_execute_gen


  

  interface fftw_plan_gen
    type(C_PTR) function fftw_plan_dft_1d(n,in,out,sign,flags) bind(C, name='fftw_plan_dft_1d')
      import
      integer(C_INT), value :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_dft_1d

    type(C_PTR) function fftw_plan_dft_2d(n0,n1,in,out,sign,flags) bind(C, name='fftw_plan_dft_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(n1,*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(n1,*), intent(inout) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_dft_2d

    type(C_PTR) function fftw_plan_dft_3d(n0,n1,n2,in,out,sign,flags) bind(C, name='fftw_plan_dft_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(n2,n1,*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(n2,n1,*), intent(inout) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_dft_3d

    type(C_PTR) function fftwf_plan_dft_1d(n,in,out,sign,flags) bind(C, name='fftwf_plan_dft_1d')
      import
      integer(C_INT), value :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_1d

    type(C_PTR) function fftwf_plan_dft_2d(n0,n1,in,out,sign,flags) bind(C, name='fftwf_plan_dft_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(n1,*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(n1,*), intent(inout) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_2d

    type(C_PTR) function fftwf_plan_dft_3d(n0,n1,n2,in,out,sign,flags) bind(C, name='fftwf_plan_dft_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(n2,n1,*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(n2,n1,*), intent(inout) :: out
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_3d

    type(C_PTR) function fftw_plan_r2r_1d(n,in,out,kind,flags) bind(C, name='fftw_plan_r2r_1d')
      import
      integer(C_INT), value :: n
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_FFTW_R2R_KIND), value :: kind
      integer(C_INT), value :: flags
    end function fftw_plan_r2r_1d

    type(C_PTR) function fftw_plan_r2r_2d(n0,n1,in,out,kind0,kind1,flags) bind(C, name='fftw_plan_r2r_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      real(C_DOUBLE), dimension(n1,*), intent(inout) :: in
      real(C_DOUBLE), dimension(n1,*), intent(inout) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftw_plan_r2r_2d

    type(C_PTR) function fftw_plan_r2r_3d(n0,n1,n2,in,out,kind0,kind1,kind2,flags) bind(C, name='fftw_plan_r2r_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      real(C_DOUBLE), dimension(n2,n1,*), intent(inout) :: in
      real(C_DOUBLE), dimension(n2,n1,*), intent(inout) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftw_plan_r2r_3d

    type(C_PTR) function fftwf_plan_r2r_1d(n,in,out,kind,flags) bind(C, name='fftwf_plan_r2r_1d')
      import
      integer(C_INT), value :: n
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_FFTW_R2R_KIND), value :: kind
      integer(C_INT), value :: flags
    end function fftwf_plan_r2r_1d

    type(C_PTR) function fftwf_plan_r2r_2d(n0,n1,in,out,kind0,kind1,flags) bind(C, name='fftwf_plan_r2r_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      real(C_FLOAT), dimension(n1,*), intent(inout) :: in
      real(C_FLOAT), dimension(n1,*), intent(inout) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftwf_plan_r2r_2d

    type(C_PTR) function fftwf_plan_r2r_3d(n0,n1,n2,in,out,kind0,kind1,kind2,flags) bind(C, name='fftwf_plan_r2r_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      real(C_FLOAT), dimension(n2,n1,*), intent(inout) :: in
      real(C_FLOAT), dimension(n2,n1,*), intent(inout) :: out
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftwf_plan_r2r_3d

    type(C_PTR) function fftw_plan_many_dft_r2c(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                         bind(C, name='fftw_plan_many_dft_r2c')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      real(C_DOUBLE), dimension(*), intent(out) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_INT), value :: flags
    end function fftw_plan_many_dft_r2c
    
    type(C_PTR) function fftw_plan_dft_r2c(rank,n,in,out,flags) bind(C, name='fftw_plan_dft_r2c')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_r2c
    
    type(C_PTR) function fftw_plan_dft_r2c_1d(n,in,out,flags) bind(C, name='fftw_plan_dft_r2c_1d')
      import
      integer(C_INT), value :: n
      real(C_DOUBLE), dimension(*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_r2c_1d
    
    type(C_PTR) function fftw_plan_dft_r2c_2d(n0,n1,in,out,flags) bind(C, name='fftw_plan_dft_r2c_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      real(C_DOUBLE), dimension(n1,*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_r2c_2d
    
    type(C_PTR) function fftw_plan_dft_r2c_3d(n0,n1,n2,in,out,flags) bind(C, name='fftw_plan_dft_r2c_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      real(C_DOUBLE), dimension(n2,n1,*), intent(out) :: in
      complex(C_DOUBLE_COMPLEX), dimension(n2,n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_r2c_3d
    
    type(C_PTR) function fftw_plan_many_dft_c2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                         bind(C, name='fftw_plan_many_dft_c2r')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_INT), value :: flags
    end function fftw_plan_many_dft_c2r
    
    type(C_PTR) function fftw_plan_dft_c2r(rank,n,in,out,flags) bind(C, name='fftw_plan_dft_c2r')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_c2r
    
    type(C_PTR) function fftw_plan_dft_c2r_1d(n,in,out,flags) bind(C, name='fftw_plan_dft_c2r_1d')
      import
      integer(C_INT), value :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_c2r_1d
    
    type(C_PTR) function fftw_plan_dft_c2r_2d(n0,n1,in,out,flags) bind(C, name='fftw_plan_dft_c2r_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(n1,*), intent(out) :: in
      real(C_DOUBLE), dimension(n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_c2r_2d
    
    type(C_PTR) function fftw_plan_dft_c2r_3d(n0,n1,n2,in,out,flags) bind(C, name='fftw_plan_dft_c2r_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(n2,n1,*), intent(out) :: in
      real(C_DOUBLE), dimension(n2,n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftw_plan_dft_c2r_3d
    
    type(C_PTR) function fftwf_plan_many_dft_r2c(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                         bind(C, name='fftwf_plan_many_dft_r2c')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      real(C_FLOAT), dimension(*), intent(out) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_INT), value :: flags
    end function fftwf_plan_many_dft_r2c
    
    type(C_PTR) function fftwf_plan_dft_r2c(rank,n,in,out,flags) bind(C, name='fftwf_plan_dft_r2c')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_r2c
    
    type(C_PTR) function fftwf_plan_dft_r2c_1d(n,in,out,flags) bind(C, name='fftwf_plan_dft_r2c_1d')
      import
      integer(C_INT), value :: n
      real(C_FLOAT), dimension(*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_r2c_1d
    
    type(C_PTR) function fftwf_plan_dft_r2c_2d(n0,n1,in,out,flags) bind(C, name='fftwf_plan_dft_r2c_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      real(C_FLOAT), dimension(n1,*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_r2c_2d
    
    type(C_PTR) function fftwf_plan_dft_r2c_3d(n0,n1,n2,in,out,flags) bind(C, name='fftwf_plan_dft_r2c_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      real(C_FLOAT), dimension(n2,n1,*), intent(out) :: in
      complex(C_FLOAT_COMPLEX), dimension(n2,n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_r2c_3d
    
    type(C_PTR) function fftwf_plan_many_dft_c2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,flags) &
                         bind(C, name='fftwf_plan_many_dft_c2r')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_INT), value :: flags
    end function fftwf_plan_many_dft_c2r
    
    type(C_PTR) function fftwf_plan_dft_c2r(rank,n,in,out,flags) bind(C, name='fftwf_plan_dft_c2r')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_c2r
    
    type(C_PTR) function fftwf_plan_dft_c2r_1d(n,in,out,flags) bind(C, name='fftwf_plan_dft_c2r_1d')
      import
      integer(C_INT), value :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_c2r_1d
    
    type(C_PTR) function fftwf_plan_dft_c2r_2d(n0,n1,in,out,flags) bind(C, name='fftwf_plan_dft_c2r_2d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(n1,*), intent(out) :: in
      real(C_FLOAT), dimension(n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_c2r_2d
    
    type(C_PTR) function fftwf_plan_dft_c2r_3d(n0,n1,n2,in,out,flags) bind(C, name='fftwf_plan_dft_c2r_3d')
      import
      integer(C_INT), value :: n0
      integer(C_INT), value :: n1
      integer(C_INT), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(n2,n1,*), intent(out) :: in
      real(C_FLOAT), dimension(n2,n1,*), intent(out) :: out
      integer(C_INT), value :: flags
    end function fftwf_plan_dft_c2r_3d
    
  end interface fftw_plan_gen

  interface
    type(C_PTR) function fftw_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags) &
                         bind(C, name='fftw_plan_many_dft')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_plan_many_dft
    
    type(C_PTR) function fftw_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags) &
                         bind(C, name='fftw_plan_many_r2r')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftw_plan_many_r2r
    
    type(C_PTR) function fftwf_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags) &
                         bind(C, name='fftwf_plan_many_dft')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_plan_many_dft
    
    type(C_PTR) function fftwf_plan_many_r2r(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,kind,flags) &
                         bind(C, name='fftwf_plan_many_r2r')
      import
      integer(C_INT), value :: rank
      integer(C_INT), dimension(*), intent(in) :: n
      integer(C_INT), value :: howmany
      real(C_FLOAT), dimension(*), intent(inout) :: in
      integer(C_INT), dimension(*), intent(in) :: inembed
      integer(C_INT), value :: istride
      integer(C_INT), value :: idist
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT), dimension(*), intent(in) :: onembed
      integer(C_INT), value :: ostride
      integer(C_INT), value :: odist
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwf_plan_many_r2r
  end interface

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
  
#ifdef MPI
  integer(C_INTPTR_T), parameter :: FFTW_MPI_DEFAULT_BLOCK = 0
  integer(C_INT), parameter :: FFTW_MPI_SCRAMBLED_IN = 134217728
  integer(C_INT), parameter :: FFTW_MPI_SCRAMBLED_OUT = 268435456
  integer(C_INT), parameter :: FFTW_MPI_TRANSPOSED_IN = 536870912
  integer(C_INT), parameter :: FFTW_MPI_TRANSPOSED_OUT = 1073741824

  type, bind(C) :: fftw_mpi_ddim
     integer(C_INTPTR_T) n, ib, ob
  end type fftw_mpi_ddim

  interface
    subroutine fftw_mpi_init() bind(C, name='fftw_mpi_init')
      import
    end subroutine fftw_mpi_init
    
    subroutine fftw_mpi_cleanup() bind(C, name='fftw_mpi_cleanup')
      import
    end subroutine fftw_mpi_cleanup
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_many_transposed(rnk,n,howmany,block0,block1,comm,local_n0,local_0_start, &
                                                                     local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_many_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_many_transposed
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_many(rnk,n,howmany,block0,comm,local_n0,local_0_start) &
                                 bind(C, name='fftw_mpi_local_size_many_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size_many
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_transposed(rnk,n,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_transposed
    
    integer(C_INTPTR_T) function fftw_mpi_local_size(rnk,n,comm,local_n0,local_0_start) bind(C, name='fftw_mpi_local_size_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_many_1d(n0,howmany,comm,sign,flags,local_ni,local_i_start,local_no, &
                                                             local_o_start) bind(C, name='fftw_mpi_local_size_many_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: howmany
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftw_mpi_local_size_many_1d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_1d(n0,comm,sign,flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='fftw_mpi_local_size_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftw_mpi_local_size_1d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_2d(n0,n1,comm,local_n0,local_0_start) &
                                 bind(C, name='fftw_mpi_local_size_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size_2d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_2d_transposed(n0,n1,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_2d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_2d_transposed
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_3d(n0,n1,n2,comm,local_n0,local_0_start) &
                                 bind(C, name='fftw_mpi_local_size_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftw_mpi_local_size_3d
    
    integer(C_INTPTR_T) function fftw_mpi_local_size_3d_transposed(n0,n1,n2,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftw_mpi_local_size_3d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftw_mpi_local_size_3d_transposed
    
    type(C_PTR) function fftw_mpi_plan_many_transpose(n0,n1,howmany,block0,block1,in,out,comm,flags) &
                         bind(C, name='fftw_mpi_plan_many_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_transpose
    
    type(C_PTR) function fftw_mpi_plan_transpose(n0,n1,in,out,comm,flags) bind(C, name='fftw_mpi_plan_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_transpose
    
    subroutine fftw_mpi_gather_wisdom(comm_) bind(C, name='fftw_mpi_gather_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftw_mpi_gather_wisdom
    
    subroutine fftw_mpi_broadcast_wisdom(comm_) bind(C, name='fftw_mpi_broadcast_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftw_mpi_broadcast_wisdom
    
    subroutine fftw_mpi_execute_dft(p,in,out) bind(C, name='fftw_mpi_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
    end subroutine fftw_mpi_execute_dft
    
    subroutine fftw_mpi_execute_dft_r2c(p,in,out) bind(C, name='fftw_mpi_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
    end subroutine fftw_mpi_execute_dft_r2c
    
    subroutine fftw_mpi_execute_dft_c2r(p,in,out) bind(C, name='fftw_mpi_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
    end subroutine fftw_mpi_execute_dft_c2r
    
    subroutine fftw_mpi_execute_r2r(p,in,out) bind(C, name='fftw_mpi_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
    end subroutine fftw_mpi_execute_r2r
    
  end interface

  type, bind(C) :: fftwf_mpi_ddim
     integer(C_INTPTR_T) n, ib, ob
  end type fftwf_mpi_ddim

  interface
    subroutine fftwf_mpi_init() bind(C, name='fftwf_mpi_init')
      import
    end subroutine fftwf_mpi_init
    
    subroutine fftwf_mpi_cleanup() bind(C, name='fftwf_mpi_cleanup')
      import
    end subroutine fftwf_mpi_cleanup
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_many_transposed(rnk,n,howmany,block0,block1,comm,local_n0,local_0_start, &
                                                                      local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_many_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_many_transposed
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_many(rnk,n,howmany,block0,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwf_mpi_local_size_many_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size_many
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_transposed(rnk,n,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_transposed_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_transposed
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size(rnk,n,comm,local_n0,local_0_start) bind(C, name='fftwf_mpi_local_size_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_many_1d(n0,howmany,comm,sign,flags,local_ni,local_i_start,local_no, &
                                                              local_o_start) bind(C, name='fftwf_mpi_local_size_many_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: howmany
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftwf_mpi_local_size_many_1d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_1d(n0,comm,sign,flags,local_ni,local_i_start,local_no,local_o_start) &
                                 bind(C, name='fftwf_mpi_local_size_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
      integer(C_INTPTR_T), intent(out) :: local_ni
      integer(C_INTPTR_T), intent(out) :: local_i_start
      integer(C_INTPTR_T), intent(out) :: local_no
      integer(C_INTPTR_T), intent(out) :: local_o_start
    end function fftwf_mpi_local_size_1d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_2d(n0,n1,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwf_mpi_local_size_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size_2d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_2d_transposed(n0,n1,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_2d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_2d_transposed
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_3d(n0,n1,n2,comm,local_n0,local_0_start) &
                                 bind(C, name='fftwf_mpi_local_size_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
    end function fftwf_mpi_local_size_3d
    
    integer(C_INTPTR_T) function fftwf_mpi_local_size_3d_transposed(n0,n1,n2,comm,local_n0,local_0_start,local_n1,local_1_start) &
                                 bind(C, name='fftwf_mpi_local_size_3d_transposed_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      integer(C_INT32_T), value :: comm
      integer(C_INTPTR_T), intent(out) :: local_n0
      integer(C_INTPTR_T), intent(out) :: local_0_start
      integer(C_INTPTR_T), intent(out) :: local_n1
      integer(C_INTPTR_T), intent(out) :: local_1_start
    end function fftwf_mpi_local_size_3d_transposed
    
    type(C_PTR) function fftwf_mpi_plan_many_transpose(n0,n1,howmany,block0,block1,in,out,comm,flags) &
                         bind(C, name='fftwf_mpi_plan_many_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block0
      integer(C_INTPTR_T), value :: block1
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_transpose
    
    type(C_PTR) function fftwf_mpi_plan_transpose(n0,n1,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_transpose_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_transpose
    

    
    subroutine fftwf_mpi_gather_wisdom(comm_) bind(C, name='fftwf_mpi_gather_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftwf_mpi_gather_wisdom
    
    subroutine fftwf_mpi_broadcast_wisdom(comm_) bind(C, name='fftwf_mpi_broadcast_wisdom_f03')
      import
      integer(C_INT32_T), value :: comm_
    end subroutine fftwf_mpi_broadcast_wisdom
    
    subroutine fftwf_mpi_execute_dft(p,in,out) bind(C, name='fftwf_mpi_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
    end subroutine fftwf_mpi_execute_dft
    
    subroutine fftwf_mpi_execute_dft_r2c(p,in,out) bind(C, name='fftwf_mpi_execute_dft_r2c')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
    end subroutine fftwf_mpi_execute_dft_r2c
    
    subroutine fftwf_mpi_execute_dft_c2r(p,in,out) bind(C, name='fftwf_mpi_execute_dft_c2r')
      import
      type(C_PTR), value :: p
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
    end subroutine fftwf_mpi_execute_dft_c2r
    
    subroutine fftwf_mpi_execute_r2r(p,in,out) bind(C, name='fftwf_mpi_execute_r2r')
      import
      type(C_PTR), value :: p
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
    end subroutine fftwf_mpi_execute_r2r
    
    
  end interface
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  interface fftw_mpi_plan_gen
  
    type(C_PTR) function fftw_mpi_plan_many_dft(rnk,n,howmany,block,tblock,in,out,comm,sign,flags) &
                         bind(C, name='fftw_mpi_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block
      integer(C_INTPTR_T), value :: tblock
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_dft
    
    type(C_PTR) function fftw_mpi_plan_dft(rnk,n,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft
    
    type(C_PTR) function fftw_mpi_plan_dft_1d(n0,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_1d
    
    type(C_PTR) function fftw_mpi_plan_dft_2d(n0,n1,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_2d
    
    type(C_PTR) function fftw_mpi_plan_dft_3d(n0,n1,n2,in,out,comm,sign,flags) bind(C, name='fftw_mpi_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_3d
    
    type(C_PTR) function fftw_mpi_plan_many_r2r(rnk,n,howmany,iblock,oblock,in,out,comm,kind,flags) &
                         bind(C, name='fftw_mpi_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_r2r
    
    type(C_PTR) function fftw_mpi_plan_r2r(rnk,n,in,out,comm,kind,flags) bind(C, name='fftw_mpi_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_r2r
    
    type(C_PTR) function fftw_mpi_plan_r2r_2d(n0,n1,in,out,comm,kind0,kind1,flags) bind(C, name='fftw_mpi_plan_r2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_r2r_2d
    
    type(C_PTR) function fftw_mpi_plan_r2r_3d(n0,n1,n2,in,out,comm,kind0,kind1,kind2,flags) bind(C, name='fftw_mpi_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_r2r_3d
    
    type(C_PTR) function fftw_mpi_plan_many_dft_r2c(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftw_mpi_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_dft_r2c
    
    type(C_PTR) function fftw_mpi_plan_dft_r2c(rnk,n,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_r2c
    
    type(C_PTR) function fftw_mpi_plan_dft_r2c_2d(n0,n1,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_r2c_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_r2c_2d
    
    type(C_PTR) function fftw_mpi_plan_dft_r2c_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_DOUBLE), dimension(*), intent(inout) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_r2c_3d
    
    type(C_PTR) function fftw_mpi_plan_many_dft_c2r(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftw_mpi_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_many_dft_c2r
    
    type(C_PTR) function fftw_mpi_plan_dft_c2r(rnk,n,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_c2r
    
    type(C_PTR) function fftw_mpi_plan_dft_c2r_2d(n0,n1,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_c2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_c2r_2d
    
    type(C_PTR) function fftw_mpi_plan_dft_c2r_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftw_mpi_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: in
      real(C_DOUBLE), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftw_mpi_plan_dft_c2r_3d
    
    
    
    type(C_PTR) function fftwf_mpi_plan_many_dft(rnk,n,howmany,block,tblock,in,out,comm,sign,flags) &
                         bind(C, name='fftwf_mpi_plan_many_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: block
      integer(C_INTPTR_T), value :: tblock
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_dft
    
    type(C_PTR) function fftwf_mpi_plan_dft(rnk,n,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft
    
    type(C_PTR) function fftwf_mpi_plan_dft_1d(n0,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_1d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_1d
    
    type(C_PTR) function fftwf_mpi_plan_dft_2d(n0,n1,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_2d
    
    type(C_PTR) function fftwf_mpi_plan_dft_3d(n0,n1,n2,in,out,comm,sign,flags) bind(C, name='fftwf_mpi_plan_dft_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: sign
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_3d
    
    type(C_PTR) function fftwf_mpi_plan_many_r2r(rnk,n,howmany,iblock,oblock,in,out,comm,kind,flags) &
                         bind(C, name='fftwf_mpi_plan_many_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_r2r
    
    type(C_PTR) function fftwf_mpi_plan_r2r(rnk,n,in,out,comm,kind,flags) bind(C, name='fftwf_mpi_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_r2r
    
    type(C_PTR) function fftwf_mpi_plan_r2r_2d(n0,n1,in,out,comm,kind0,kind1,flags) bind(C, name='fftwf_mpi_plan_r2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_r2r_2d
    
    type(C_PTR) function fftwf_mpi_plan_r2r_3d(n0,n1,n2,in,out,comm,kind0,kind1,kind2,flags) &
                         bind(C, name='fftwf_mpi_plan_r2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_FLOAT), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_FFTW_R2R_KIND), value :: kind0
      integer(C_FFTW_R2R_KIND), value :: kind1
      integer(C_FFTW_R2R_KIND), value :: kind2
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_r2r_3d
    
    type(C_PTR) function fftwf_mpi_plan_many_dft_r2c(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftwf_mpi_plan_many_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_dft_r2c
    
    type(C_PTR) function fftwf_mpi_plan_dft_r2c(rnk,n,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_r2c_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_r2c
    
    type(C_PTR) function fftwf_mpi_plan_dft_r2c_2d(n0,n1,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_r2c_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_r2c_2d
    
    type(C_PTR) function fftwf_mpi_plan_dft_r2c_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_r2c_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      real(C_FLOAT), dimension(*), intent(inout) :: in
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_r2c_3d
    
    type(C_PTR) function fftwf_mpi_plan_many_dft_c2r(rnk,n,howmany,iblock,oblock,in,out,comm,flags) &
                         bind(C, name='fftwf_mpi_plan_many_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      integer(C_INTPTR_T), value :: howmany
      integer(C_INTPTR_T), value :: iblock
      integer(C_INTPTR_T), value :: oblock
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_many_dft_c2r
    
    type(C_PTR) function fftwf_mpi_plan_dft_c2r(rnk,n,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_c2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: n
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_c2r
    
    type(C_PTR) function fftwf_mpi_plan_dft_c2r_2d(n0,n1,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_c2r_2d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_c2r_2d
    
    type(C_PTR) function fftwf_mpi_plan_dft_c2r_3d(n0,n1,n2,in,out,comm,flags) bind(C, name='fftwf_mpi_plan_dft_c2r_3d_f03')
      import
      integer(C_INTPTR_T), value :: n0
      integer(C_INTPTR_T), value :: n1
      integer(C_INTPTR_T), value :: n2
      complex(C_FLOAT_COMPLEX), dimension(*), intent(inout) :: in
      real(C_FLOAT), dimension(*), intent(inout) :: out
      integer(C_INT32_T), value :: comm
      integer(C_INT), value :: flags
    end function fftwf_mpi_plan_dft_c2r_3d
  end interface fftw_mpi_plan_gen
#endif
end module fftw3