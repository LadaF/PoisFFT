module PFFT
#ifdef MPI
  use fftw3
  
#include "pfft.f03"
! 
! ! Generated automatically.  DO NOT EDIT! :) :(
! 
! ! integers
! 
! ! unsigned
!   integer(C_INT), parameter :: PFFT_TRANSPOSED_NONE = 0
!   integer(C_INT), parameter :: PFFT_MEASURE = 0
!   integer(C_INT), parameter :: PFFT_NO_TUNE = 0
!   integer(C_INT), parameter :: FPFFT_DEFAULT_BLOCKS = -1
!   integer(C_INT), parameter :: FPFFT_NO_GCELLS = -1
!   integer(C_INT), parameter :: PFFT_INT = 1
!   integer(C_INT), parameter :: PFFT_PTRDIFF_T = 2
!   integer(C_INT), parameter :: PFFT_FLOAT = 3
!   integer(C_INT), parameter :: PFFT_DOUBLE = 4
!   integer(C_INT), parameter :: PFFT_UNSIGNED = 5
!   integer(C_INT), parameter :: PFFT_GC_NONTRANSPOSED = 0
! 
! ! shifted unsigned
!   integer(C_INT), parameter :: PFFT_TRANSPOSED_IN = 1
!   integer(C_INT), parameter :: PFFT_TRANSPOSED_OUT = 2
!   integer(C_INT), parameter :: PFFT_SHIFTED_IN = 4
!   integer(C_INT), parameter :: PFFT_SHIFTED_OUT = 8
!   integer(C_INT), parameter :: PFFT_ESTIMATE = 16
!   integer(C_INT), parameter :: PFFT_PATIENT = 32
!   integer(C_INT), parameter :: PFFT_EXHAUSTIVE = 64
!   integer(C_INT), parameter :: PFFT_TUNE = 128
!   integer(C_INT), parameter :: PFFT_PRESERVE_INPUT = 256
!   integer(C_INT), parameter :: PFFT_DESTROY_INPUT = 512
!   integer(C_INT), parameter :: PFFT_BUFFERED_INPLACE = 1024
!   integer(C_INT), parameter :: PFFT_GC_TRANSPOSED = 1
!   integer(C_INT), parameter :: PFFT_GC_SENDRECV = 2
!   integer(C_INT), parameter :: PFFT_GC_RMA = 4
! 
! ! redirections
!   integer(C_INT), parameter :: PFFT_R2HC = FFTW_R2HC
!   integer(C_INT), parameter :: PFFT_HC2R = FFTW_HC2R
!   integer(C_INT), parameter :: PFFT_DHT = FFTW_DHT
!   integer(C_INT), parameter :: PFFT_REDFT00 = FFTW_REDFT00
!   integer(C_INT), parameter :: PFFT_REDFT01 = FFTW_REDFT01
!   integer(C_INT), parameter :: PFFT_REDFT10 = FFTW_REDFT10
!   integer(C_INT), parameter :: PFFT_REDFT11 = FFTW_REDFT11
!   integer(C_INT), parameter :: PFFT_RODFT00 = FFTW_RODFT00
!   integer(C_INT), parameter :: PFFT_RODFT01 = FFTW_RODFT01
!   integer(C_INT), parameter :: PFFT_RODFT10 = FFTW_RODFT10
!   integer(C_INT), parameter :: PFFT_RODFT11 = FFTW_RODFT11
!   integer(C_INT), parameter :: PFFT_FORWARD = FFTW_FORWARD
!   integer(C_INT), parameter :: PFFT_BACKWARD = FFTW_BACKWARD
!   integer(C_INT), parameter :: PFFT_DEFAULT_BLOCK = FFTW_MPI_DEFAULT_BLOCK
!   
!   interface pfft_plan_gen
!   
!     type(C_PTR) function pfft_plan_dft_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_dft_3d
!     
!    type(C_PTR) function pfft_plan_dft_r2c_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_r2c_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       real(C_DOUBLE), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_dft_r2c_3d
!     
!     type(C_PTR) function pfft_plan_dft_c2r_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_c2r_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       real(C_DOUBLE), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_dft_c2r_3d
!     
!     type(C_PTR) function pfft_plan_r2r_3d(Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfft_plan_r2r_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       real(C_DOUBLE), dimension(*), intent(out) :: in
!       real(C_DOUBLE), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_r2r_3d
!     
!     type(C_PTR) function pfft_plan_dft(rnk,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_f03')
!       import
!       integer(C_INT), value :: rnk
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_dft
!     
!     type(C_PTR) function pfft_plan_dft_r2c(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_r2c_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       real(C_DOUBLE), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_dft_r2c
!     
!     type(C_PTR) function pfft_plan_dft_c2r(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfft_plan_dft_c2r_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       real(C_DOUBLE), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_dft_c2r
!     
!     type(C_PTR) function pfft_plan_many_dft(rnk,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
!                          bind(C, name='pfft_plan_many_dft_f03')
!       import
!       integer(C_INT), value :: rnk
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_many_dft
!     
!     type(C_PTR) function pfft_plan_many_dft_r2c(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
!                          bind(C, name='pfft_plan_many_dft_r2c_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       real(C_DOUBLE), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_many_dft_r2c
!     
!     type(C_PTR) function pfft_plan_many_dft_c2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
!                          bind(C, name='pfft_plan_many_dft_c2r_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       real(C_DOUBLE), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_many_dft_c2r
!     
!     type(C_PTR) function pfft_plan_many_r2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,kinds,pfft_flags) &
!                          bind(C, name='pfft_plan_many_r2r_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       real(C_DOUBLE), dimension(*), intent(out) :: in
!       real(C_DOUBLE), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_many_r2r
!     
!     type(C_PTR) function pfft_plan_many_dft_skipped(rnk,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
!                                                     pfft_flags) bind(C, name='pfft_plan_many_dft_skipped_f03')
!       import
!       integer(C_INT), value :: rnk
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       integer(C_INT), dimension(*), intent(in) :: skip_trafos
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfft_plan_many_dft_skipped
!     
!     
!     
!     
!     type(C_PTR) function pfftf_plan_dft_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_dft_3d
!     
!     type(C_PTR) function pfftf_plan_dft_r2c_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_r2c_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       real(C_FLOAT), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_dft_r2c_3d
!     
!     type(C_PTR) function pfftf_plan_dft_c2r_3d(Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_c2r_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       real(C_FLOAT), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_dft_c2r_3d
!     
!     type(C_PTR) function pfftf_plan_r2r_3d(Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfftf_plan_r2r_3d_f03')
!       import
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       real(C_FLOAT), dimension(*), intent(out) :: in
!       real(C_FLOAT), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_r2r_3d
!     
!     type(C_PTR) function pfftf_plan_dft(rnk,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_f03')
!       import
!       integer(C_INT), value :: rnk
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_dft
!     
!     type(C_PTR) function pfftf_plan_dft_r2c(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_r2c_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       real(C_FLOAT), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_dft_r2c
!     
!     type(C_PTR) function pfftf_plan_dft_c2r(rnk_n,Nos,in,out,comm_cart,sign,pfft_flags) bind(C, name='pfftf_plan_dft_c2r_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       real(C_FLOAT), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_dft_c2r
!     
!     type(C_PTR) function pfftf_plan_many_dft(rnk,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
!                          bind(C, name='pfftf_plan_many_dft_f03')
!       import
!       integer(C_INT), value :: rnk
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_many_dft
!     
!     type(C_PTR) function pfftf_plan_many_dft_r2c(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
!                          bind(C, name='pfftf_plan_many_dft_r2c_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       real(C_FLOAT), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_many_dft_r2c
!     
!     type(C_PTR) function pfftf_plan_many_dft_c2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,sign,pfft_flags) &
!                          bind(C, name='pfftf_plan_many_dft_c2r_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       real(C_FLOAT), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_many_dft_c2r
!     
!     type(C_PTR) function pfftf_plan_many_r2r(rnk_n,Nos,ni,no,howmany,iblock,oblock,in,out,comm_cart,kinds,pfft_flags) &
!                          bind(C, name='pfftf_plan_many_r2r_f03')
!       import
!       integer(C_INT), value :: rnk_n
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       real(C_FLOAT), dimension(*), intent(out) :: in
!       real(C_FLOAT), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_many_r2r
!     
!     type(C_PTR) function pfftf_plan_many_dft_skipped(rnk,Nos,ni,no,howmany,iblock,oblock,skip_trafos,in,out,comm_cart,sign, &
!                                                      pfft_flags) bind(C, name='pfftf_plan_many_dft_skipped_f03')
!       import
!       integer(C_INT), value :: rnk
!       integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
!       integer(C_INTPTR_T), dimension(*), intent(in) :: ni
!       integer(C_INTPTR_T), dimension(*), intent(in) :: no
!       integer(C_INTPTR_T), value :: howmany
!       integer(C_INTPTR_T), dimension(*), intent(in) :: iblock
!       integer(C_INTPTR_T), dimension(*), intent(in) :: oblock
!       integer(C_INT), dimension(*), intent(in) :: skip_trafos
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: in
!       complex(C_FLOAT_COMPLEX), dimension(*), intent(out) :: out
!       integer(C_INT32_T), value :: comm_cart
!       integer(C_INT), value :: sign
!       integer(C_INT), value :: pfft_flags
!     end function pfftf_plan_many_dft_skipped
    

!   end interface
#endif
end module PFFT

