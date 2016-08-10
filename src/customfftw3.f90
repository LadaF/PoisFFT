module fftw3
use iso_c_binding

#include <fftw3.f03>
 
  interface fftw_execute_gen
    procedure fftw_execute_r2r
    procedure fftwf_execute_r2r
    procedure fftw_execute_dft
    procedure fftwf_execute_dft
    procedure fftw_execute_dft_r2c
    procedure fftwf_execute_dft_r2c
    procedure fftw_execute_dft_c2r
    procedure fftwf_execute_dft_c2r
  end interface fftw_execute_gen

  interface fftw_plan_gen
    procedure fftw_plan_dft_1d
    procedure fftw_plan_dft_2d
    procedure fftw_plan_dft_3d
    procedure fftwf_plan_dft_1d
    procedure fftwf_plan_dft_2d
    procedure fftwf_plan_dft_3d
    procedure fftw_plan_r2r_1d
    procedure fftw_plan_r2r_2d
    procedure fftw_plan_r2r_3d
    procedure fftwf_plan_r2r_1d
    procedure fftwf_plan_r2r_2d
    procedure fftwf_plan_r2r_3d
    procedure fftw_plan_dft_r2c_1d
    procedure fftw_plan_dft_r2c_2d
    procedure fftw_plan_dft_r2c_3d
    procedure fftw_plan_dft_c2r_1d
    procedure fftw_plan_dft_c2r_2d
    procedure fftw_plan_dft_c2r_3d
    procedure fftwf_plan_dft_r2c_1d
    procedure fftwf_plan_dft_r2c_2d
    procedure fftwf_plan_dft_r2c_3d
    procedure fftwf_plan_dft_c2r_1d
    procedure fftwf_plan_dft_c2r_2d
    procedure fftwf_plan_dft_c2r_3d
  end interface fftw_plan_gen


  
#ifdef MPI

#include <fftw3-mpi.f03>
  
  interface fftw_mpi_plan_gen
  
    procedure fftw_mpi_plan_many_dft
    
    procedure fftw_mpi_plan_dft
    
    procedure fftw_mpi_plan_dft_1d
    
    procedure fftw_mpi_plan_dft_2d
    
    procedure fftw_mpi_plan_dft_3d
    
    procedure fftw_mpi_plan_many_r2r
    
    procedure fftw_mpi_plan_r2r
    
    procedure fftw_mpi_plan_r2r_2d
    
    procedure fftw_mpi_plan_r2r_3d
    
    procedure fftw_mpi_plan_many_dft_r2c
    
    procedure fftw_mpi_plan_dft_r2c
    
    procedure fftw_mpi_plan_dft_r2c_2d
    
    procedure fftw_mpi_plan_dft_r2c_3d
    
    procedure fftw_mpi_plan_many_dft_c2r
    
    procedure fftw_mpi_plan_dft_c2r
    
    procedure fftw_mpi_plan_dft_c2r_2d
    
    procedure fftw_mpi_plan_dft_c2r_3d
    
    
    
    procedure fftwf_mpi_plan_many_dft
    
    procedure fftwf_mpi_plan_dft
    
    procedure fftwf_mpi_plan_dft_1d
    
    procedure fftwf_mpi_plan_dft_2d
    
    procedure fftwf_mpi_plan_dft_3d
    
    procedure fftwf_mpi_plan_many_r2r
    
    procedure fftwf_mpi_plan_r2r
    
    procedure fftwf_mpi_plan_r2r_2d
    
    procedure fftwf_mpi_plan_r2r_3d
    
    procedure fftwf_mpi_plan_many_dft_r2c
    
    procedure fftwf_mpi_plan_dft_r2c
    
    procedure fftwf_mpi_plan_dft_r2c_2d
    
    procedure fftwf_mpi_plan_dft_r2c_3d
    
    procedure fftwf_mpi_plan_many_dft_c2r
    
    procedure fftwf_mpi_plan_dft_c2r
    
    procedure fftwf_mpi_plan_dft_c2r_2d
    
    procedure fftwf_mpi_plan_dft_c2r_3d
  end interface fftw_mpi_plan_gen
#endif
end module fftw3