module PFFT
#ifdef MPI
  use fftw3
  
#include "pfft.f03"

!NOTE: only for older versions of PFFT
#ifdef MISSING_PFFT_PLAN_WITH_NTHREADS
  interface
    subroutine pfft_plan_with_nthreads(nthreads) bind(C, name="pfft_plan_with_nthreads")
      import
      integer(C_INT), value :: nthreads
    end subroutine
  end interface
#endif

!NOTE: only for older versions of PFFT
#ifdef MISSING_PFFT_R2R
  interface
    type(C_PTR) function pfft_plan_r2r(rnk,Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfft_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_DOUBLE), dimension(*), intent(out) :: in
      real(C_DOUBLE), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfft_plan_r2r
    type(C_PTR) function pfftf_plan_r2r(rnk,Nos,in,out,comm_cart,kinds,pfft_flags) bind(C, name='pfftf_plan_r2r_f03')
      import
      integer(C_INT), value :: rnk
      integer(C_INTPTR_T), dimension(*), intent(in) :: Nos
      real(C_FLOAT), dimension(*), intent(out) :: in
      real(C_FLOAT), dimension(*), intent(out) :: out
      integer(C_INT32_T), value :: comm_cart
      integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kinds
      integer(C_INT), value :: pfft_flags
    end function pfftf_plan_r2r
  end interface
#endif
  
#endif
end module PFFT

