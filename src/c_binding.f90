module PoisFFT_C_binding
  use iso_c_binding
  use poisfft
  
  implicit none
  
#ifdef MPI
  interface
    integer function MPI_Comm_c2f(c_handle) bind(C, name="f_MPI_Comm_c2f")
      use iso_c_binding
      type(c_ptr), value :: c_handle
    end function
  end interface
#endif

contains

#define rp c_double
  subroutine poisfft_solver1d_new(D, nxyz, Lxyz, BCs, approximation, &
                                  gnxyz, offs, comm_ptr, nthreads) &
    bind(C, name="poisfft_solver1d_new")
#define dims 1
#define solver PoisFFT_Solver1D_DP
#include "c_new-inc.f90"
#undef solver
#undef dims
  end subroutine
  
  subroutine poisfft_solver2d_new(D, nxyz, Lxyz, BCs, approximation, &
                                  gnxyz, offs, comm_ptr, nthreads) &
    bind(C, name="poisfft_solver2d_new")
#define dims 2
#define solver PoisFFT_Solver2D_DP
#include "c_new-inc.f90"
#undef solver
#undef dims
  end subroutine
  
  subroutine poisfft_solver3d_new(D, nxyz, Lxyz, BCs, approximation, &
                                  gnxyz, offs, comm_ptr, nthreads) &
    bind(C, name="poisfft_solver3d_new")
#define dims 3
#define solver PoisFFT_Solver3D_DP
#include "c_new-inc.f90"
#undef solver
#undef dims
  end subroutine
#undef rp
  
#define rp c_float
  subroutine poisfft_solver1d_f_new(D, nxyz, Lxyz, BCs, approximation, &
                                    gnxyz, offs, comm_ptr, nthreads) &
    bind(C, name="poisfft_solver1d_f_new")
#define dims 1
#define solver PoisFFT_Solver1D_SP
#include "c_new-inc.f90"
#undef solver
#undef dims
  end subroutine
  
  subroutine poisfft_solver2d_f_new(D, nxyz, Lxyz, BCs, approximation, &
                                    gnxyz, offs, comm_ptr, nthreads) &
    bind(C, name="poisfft_solver2d_f_new")
#define dims 2
#define solver PoisFFT_Solver2D_SP
#include "c_new-inc.f90"
#undef solver
#undef dims
  end subroutine
  
  subroutine poisfft_solver3d_f_new(D, nxyz, Lxyz, BCs, approximation, &
                                    gnxyz, offs, comm_ptr, nthreads) &
    bind(C, name="poisfft_solver3d_f_new")
#define dims 3
#define solver PoisFFT_Solver3D_SP
#include "c_new-inc.f90"
#undef solver
#undef dims
  end subroutine
#undef rp



#define rp c_double    
  subroutine poisfft_solver1d_execute(D, Phi, RHS, ngPhi, ngRHS) &
    bind(C, name="poisfft_solver1d_execute")
    type(PoisFFT_Solver1D_DP), pointer :: f_D
#define dims 1
#include "c_execute-inc.f90"
#undef dims
  end subroutine

  subroutine poisfft_solver2d_execute(D, Phi, RHS, ngPhi, ngRHS) &
    bind(C, name="poisfft_solver2d_execute")
    type(PoisFFT_Solver2D_DP), pointer :: f_D
#define dims 2
#include "c_execute-inc.f90"
#undef dims
  end subroutine

  subroutine poisfft_solver3d_execute(D, Phi, RHS, ngPhi, ngRHS) &
    bind(C, name="poisfft_solver3d_execute")
    type(PoisFFT_Solver3D_DP), pointer :: f_D
#define dims 3
#include "c_execute-inc.f90"
#undef dims
  end subroutine
#undef dims
#undef rp

#define rp c_float    
  subroutine poisfft_solver1d_f_execute(D, Phi, RHS, ngPhi, ngRHS) &
    bind(C, name="poisfft_solver1d_f_execute")
    type(PoisFFT_Solver1D_SP), pointer :: f_D
#define dims 1
#include "c_execute-inc.f90"
#undef dims
  end subroutine

  subroutine poisfft_solver2d_f_execute(D, Phi, RHS, ngPhi, ngRHS) &
    bind(C, name="poisfft_solver2d_f_execute")
    type(PoisFFT_Solver2D_SP), pointer :: f_D
#define dims 2
#include "c_execute-inc.f90"
#undef dims
  end subroutine

  subroutine poisfft_solver3d_f_execute(D, Phi, RHS, ngPhi, ngRHS) &
    bind(C, name="poisfft_solver3d_f_execute")
    type(PoisFFT_Solver3D_SP), pointer :: f_D
#define dims 3
#include "c_execute-inc.f90"
#undef dims
  end subroutine
#undef dims
#undef rp


  subroutine poisfft_solver1d_finalize(D) &
    bind(C, name="poisfft_solver1d_finalize")
    type(PoisFFT_Solver1D_DP), pointer :: f_D
#include "c_finalize-inc.f90"
  end subroutine

  subroutine poisfft_solver2d_finalize(D) &
    bind(C, name="poisfft_solver2d_finalize")
    type(PoisFFT_Solver2D_DP), pointer :: f_D
#include "c_finalize-inc.f90"
  end subroutine

  subroutine poisfft_solver3d_finalize(D) &
    bind(C, name="poisfft_solver3d_finalize")
    type(PoisFFT_Solver3D_DP), pointer :: f_D
#include "c_finalize-inc.f90"
  end subroutine

  subroutine poisfft_solver1d_f_finalize(D) &
    bind(C, name="poisfft_solver1d_f_finalize")
    type(PoisFFT_Solver1D_SP), pointer :: f_D
#include "c_finalize-inc.f90"
  end subroutine

  subroutine poisfft_solver2d_f_finalize(D) &
    bind(C, name="poisfft_solver2d_f_finalize")
    type(PoisFFT_Solver2D_SP), pointer :: f_D
#include "c_finalize-inc.f90"
  end subroutine

  subroutine poisfft_solver3d_f_finalize(D) &
    bind(C, name="poisfft_solver3d_f_finalize")
    type(PoisFFT_Solver3D_SP), pointer :: f_D
#include "c_finalize-inc.f90"
  end subroutine

end module PoisFFT_C_binding
