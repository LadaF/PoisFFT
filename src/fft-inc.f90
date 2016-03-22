#if (PREC==2)
#define RP DRP
#define CP DCP
#define fftw_execute_complex fftw_execute_dft
#define fftw_mpi_execute_complex fftw_mpi_execute_dft
#define fftw_execute_real fftw_execute_r2r
#define pfft_execute_gen pfft_execute
#else
#define RP SRP
#define CP SCP
#define fftw_execute_complex fftwf_execute_dft
#define fftw_mpi_execute_complex fftwf_mpi_execute_dft
#define fftw_execute_real fftwf_execute_r2r
#define pfft_execute_gen pfftf_execute
#endif
  use iso_c_binding
  use PoisFFT_Precisions
  use fftw3
#ifdef MPI
  use pfft
#endif
  implicit none



  integer, parameter :: FFT_Complex = 0
  integer, parameter :: FFT_RealEven00  = FFTW_REDFT00
  integer, parameter :: FFT_RealEven01  = FFTW_REDFT01
  integer, parameter :: FFT_RealEven10  = FFTW_REDFT10
  integer, parameter :: FFT_RealEven11  = FFTW_REDFT11
  integer, parameter :: FFT_RealOdd00   = FFTW_RODFT00
  integer, parameter :: FFT_RealOdd01   = FFTW_RODFT01
  integer, parameter :: FFT_RealOdd10   = FFTW_RODFT10
  integer, parameter :: FFT_RealOdd11   = FFTW_RODFT11

  integer, parameter :: FFT_DISTRIBUTED_FFTW = 1
  integer, parameter :: FFT_DISTRIBUTED_PFFT = 2

  type PoisFFT_Plan1D
    type(c_ptr)         :: planptr = c_null_ptr
    logical             :: planowner = .false.
    logical             :: distributed = .false.
    integer(c_int)      :: dir
  end type PoisFFT_Plan1D


  type PoisFFT_Plan2D
    type(c_ptr)         :: planptr = c_null_ptr
    logical             :: planowner = .false.
    logical             :: distributed = .false.
    integer             :: method = FFT_DISTRIBUTED_PFFT
    integer(c_int)      :: dir
  end type PoisFFT_Plan2D


  type PoisFFT_Plan3D
    type(c_ptr)         :: planptr = c_null_ptr
    logical             :: planowner = .false.
    logical             :: distributed = .false.
    integer(c_int)      :: dir
  end type PoisFFT_Plan3D


  type mpi_vars_1d
    !MPI communicator for the exchange
    integer :: comm = -1
    !our rank in comm
    integer :: rank
    !number of processes in comm
    integer :: np
    ! s for dimensions of send buffers, r for receive buffers
    ! nx, nz are the dimensions of individual parts which will be sent to or received from other processes in comm
    ! displs are the displacements (offsets) of individual pieces in the whole 1D buffer
    ! counts are the the number of elements in each piece
    ! sumrnzs(i) is the sum of rnzs in pieces 1..i-1
    integer,dimension(:),allocatable :: snxs, snzs, &
                                        rnxs, rnzs, &
                                        sdispls, scounts, &
                                        rdispls, rcounts, &
                                        sumrnzs
    !contiguous 3D buffer for locally transposed RHS
    real(RP), allocatable :: tmp1(:,:,:)
    !1D buffer for globally transposed blocks of tmp1 from different processes after MPI_Alltoallv
    real(RP), allocatable :: tmp2(:)
    !3D buffer for contiguous 1D rows on which 1D FFT can be performed locally,
    ! constructed by local reordering of tmp2
    real(RP), allocatable :: rwork(:,:,:)
  end type

  type PoisFFT_Solver1D
    real(RP) :: Lx
    integer(c_int) :: nxyz(1)
    integer(c_int) :: nx
    integer(c_int) :: gnx
    integer(c_int) :: offx = 0 !offset from global index
    integer(c_size_t) :: cnt
    integer(c_size_t) :: gcnt
    real(RP) :: norm_factor
    integer(c_int), dimension(2) :: BCs
    real(RP), allocatable, dimension(:) :: denomx
    integer :: approximation = 0
    
    type(PoisFFT_Plan1D) :: forward, backward
    complex(CP), dimension(:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: cwork => null()
    real(RP), dimension(:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: rwork => null()
    logical :: mpi_transpose_needed = .false.
    integer :: nthreads = 1
    type(mpi_vars_1D) :: mpi
  end type PoisFFT_Solver1D

  type mpi_vars_2d
    integer :: rank
    integer :: np
    integer :: comm = -1
    integer :: comm_dim = 2
  end type

  type PoisFFT_Solver2D
    real(RP) :: Lx, Ly
    integer(c_int) :: nxyz(2)
    integer(c_int) :: nx, ny
    integer(c_int) :: gnx, gny
    integer(c_int) :: offx = 0, offy = 0 !offsets from global indexes
    integer(c_size_t) :: cnt
    integer(c_size_t) :: gcnt
    real(RP) :: norm_factor
    integer(c_int), dimension(4) :: BCs
    real(RP), allocatable, dimension(:) :: denomx, denomy
    integer :: approximation = 0
    
    type(PoisFFT_Plan2D) :: forward, backward
    complex(CP), dimension(:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: cwork => null()
    real(RP), dimension(:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: rwork => null()
    integer :: nthreads = 1
    type(mpi_vars_2D) :: mpi
  end type PoisFFT_Solver2D
  
  type mpi_vars_3D
    integer :: comm = -1
  end type

  type PoisFFT_Solver3D
    real(RP) :: Lx, Ly, Lz
    integer(c_int) :: nxyz(3)
    integer(c_int) :: nx, ny, nz
    integer(c_int) :: gnx, gny, gnz
    integer(c_int) :: offx = 0, offy = 0, offz = 0 !offsets from global indexes
    integer(c_size_t) :: cnt
    integer(c_size_t) :: gcnt
    real(RP) :: norm_factor
    integer(c_int), dimension(6) :: BCs
    real(RP), allocatable, dimension(:) :: denomx, denomy, denomz
    integer :: approximation = 0
    
    type(PoisFFT_Plan3D) :: forward, backward
    complex(CP), dimension(:,:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: cwork => null()
    real(RP), dimension(:,:,:), &
#ifndef NO_CONTIGUOUS
    contiguous, &
#endif
                       pointer :: rwork => null()
    integer :: nthreads = 1
    !will be used in splitting for some boundary conditions
    type(PoisFFT_Solver1D),dimension(:),allocatable :: Solvers1D
    type(PoisFFT_Solver2D),dimension(:),allocatable :: Solvers2D
    type(mpi_vars_3D) :: mpi
 end type PoisFFT_Solver3D

  interface deallocate_fftw
    module procedure deallocate_fftw_1D_complex
    module procedure deallocate_fftw_1D_real
    module procedure deallocate_fftw_2D_complex
    module procedure deallocate_fftw_2D_real
    module procedure deallocate_fftw_3D_complex
    module procedure deallocate_fftw_3D_real
  end interface

  interface allocate_fftw_real
    module procedure allocate_fftw_1D_real
    module procedure allocate_fftw_2D_real
    module procedure allocate_fftw_3D_real
  end interface

  interface allocate_fftw_complex
    module procedure allocate_fftw_1D_complex
    module procedure allocate_fftw_2D_complex
    module procedure allocate_fftw_3D_complex
  end interface

  interface Execute
    module procedure PoisFFT_Plan1D_Execute_Real
    module procedure PoisFFT_Plan1D_Execute_Complex
    module procedure PoisFFT_Plan2D_Execute_Real
    module procedure PoisFFT_Plan2D_Execute_Complex
    module procedure PoisFFT_Plan3D_Execute_Real
    module procedure PoisFFT_Plan3D_Execute_Complex
  end interface
  
#ifdef MPI
  interface Execute_MPI
    module procedure PoisFFT_Plan1D_Execute_MPI
    module procedure PoisFFT_Plan2D_Execute_MPI
    module procedure PoisFFT_Plan3D_Execute_MPI
 end interface
#endif

  interface Finalize
    module procedure PoisFFT_Plan1D_Finalize
    module procedure PoisFFT_Plan2D_Finalize
    module procedure PoisFFT_Plan3D_Finalize
  end interface

  interface PoisFFT_Plan1D
    module procedure PoisFFT_Plan1D_New
  end interface
  interface PoisFFT_Plan2D
    module procedure PoisFFT_Plan2D_New
  end interface
  interface PoisFFT_Plan3D
    module procedure PoisFFT_Plan3D_New
  end interface


  contains


    function PoisFFT_Plan1D_New(D, plantypes, distributed) result(plan)
#define dimensions 1
#include "plan_new-inc.f90"
#undef dimensions
    end function

    function PoisFFT_Plan2D_New(D, plantypes, distributed) result(plan)
#define dimensions 2
#include "plan_new-inc.f90"
#undef dimensions
    end function

    function PoisFFT_Plan3D_New(D, plantypes, distributed) result(plan)
#define dimensions 3
#include "plan_new-inc.f90"
#undef dimensions
    end function






    subroutine allocate_fftw_1D_complex(D)
#define dimensions 1
#define realcomplex 2

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
    end subroutine

    subroutine allocate_fftw_1D_real(D)
#define dimensions 1
#define realcomplex 1

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
    end subroutine



    subroutine allocate_fftw_2D_complex(D)
#define dimensions 2
#define realcomplex 2

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
    end subroutine

    subroutine allocate_fftw_2D_real(D)
#define dimensions 2
#define realcomplex 1

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
    end subroutine

    
    

    subroutine allocate_fftw_3D_complex(D)
#define dimensions 3
#define realcomplex 2

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
    end subroutine

    subroutine allocate_fftw_3D_real(D)
#define dimensions 3
#define realcomplex 1

#include "allocate_fftw-inc.f90"

#undef dimensions
#undef realcomplex
    end subroutine




    subroutine PoisFFT_Plan1D_Execute_Complex(plan, data)
      type(PoisFFT_Plan1D), intent(in) :: plan
      complex(CP), dimension(:) &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      call fftw_execute_complex(plan%planptr, data, data)

    end subroutine

    subroutine PoisFFT_Plan1D_Execute_Real(plan, data)
      type(PoisFFT_Plan1D), intent(in) :: plan
      real(RP), dimension(:) &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      call fftw_execute_real(plan%planptr, data, data)

    end subroutine
    

    subroutine PoisFFT_Plan2D_Execute_Complex(plan, data)
      type(PoisFFT_Plan2D), intent(in) :: plan
      complex(CP), dimension(:,:) &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      call fftw_execute_complex(plan%planptr, data, data)

    end subroutine


    subroutine PoisFFT_Plan2D_Execute_Real(plan, data)
      type(PoisFFT_Plan2D), intent(in) :: plan
      real(RP), dimension(:,:) &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      call fftw_execute_real(plan%planptr, data, data)

    end subroutine
    

    subroutine PoisFFT_Plan3D_Execute_Complex(plan, data)
      type(PoisFFT_Plan3D), intent(in) :: plan
      complex(CP), dimension(:,:,:) &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      call fftw_execute_complex(plan%planptr, data, data)

    end subroutine

    subroutine PoisFFT_Plan3D_Execute_Real(plan, data)
      type(PoisFFT_Plan3D), intent(in) :: plan
      real(RP), dimension(:,:,:) &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      call fftw_execute_real(plan%planptr, data, data)

    end subroutine



  
#ifdef MPI
    subroutine PoisFFT_Plan1D_Execute_MPI(plan)
      type(PoisFFT_Plan1D), intent(in) :: plan

        call pfft_execute_gen(plan%planptr)

    end subroutine

    subroutine PoisFFT_Plan2D_Execute_MPI(plan, data)
      type(PoisFFT_Plan2D), intent(in) :: plan
      complex(CP), dimension(:,:), optional &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data

      if (plan%method==FFT_DISTRIBUTED_FFTW) then
        if (present(data)) then
          call fftw_mpi_execute_complex(plan%planptr, data, data)
        else
          stop "Error, missing `data` argument in a call to Execute_MPI with FFTW."
        end if
      else
        call pfft_execute_gen(plan%planptr)
      end if

    end subroutine

    subroutine PoisFFT_Plan3D_Execute_MPI(plan)
      type(PoisFFT_Plan3D), intent(in) :: plan

        call pfft_execute_gen(plan%planptr)

    end subroutine

#endif
    
    
    
    

    subroutine deallocate_fftw_1D_complex(data)
      complex(CP), dimension(:), pointer &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data
      type(c_ptr) :: p

      p = c_loc(data(1))
      call fftw_free(p)
      nullify(data)
    end subroutine deallocate_fftw_1D_complex

    subroutine deallocate_fftw_1D_real(data)
      real(RP), dimension(:), pointer &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data
      type(c_ptr) :: p

      p = c_loc(data(1))
      call fftw_free(p)
      nullify(data)
    end subroutine deallocate_fftw_1D_real

    subroutine deallocate_fftw_2D_complex(data)
      complex(CP), dimension(:,:), pointer &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine deallocate_fftw_2D_complex

    subroutine deallocate_fftw_2D_real(data)
      real(RP), dimension(:,:), pointer &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine deallocate_fftw_2D_real

    subroutine deallocate_fftw_3D_complex(data)
      complex(CP), dimension(:,:,:), pointer &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine deallocate_fftw_3D_complex


    subroutine deallocate_fftw_3D_real(data)
      real(RP), dimension(:,:,:), pointer &
#ifndef NO_CONTIGUOUS
                               , contiguous &
#endif
                                           :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine deallocate_fftw_3D_real





    subroutine PoisFFT_Plan1D_Finalize(plan)
      type(PoisFFT_Plan1D) :: plan

      if (c_associated(plan%planptr).and.plan%planowner)  &
        call fftw_destroy_plan(plan%planptr)
      plan%planptr = c_null_ptr
    end subroutine PoisFFT_Plan1D_Finalize

    subroutine PoisFFT_Plan2D_Finalize(plan)
      type(PoisFFT_Plan2D) :: plan

      if (c_associated(plan%planptr).and.plan%planowner) then
#ifdef MPI
        if (plan%distributed.and.plan%method==FFT_DISTRIBUTED_PFFT) then
          call pfft_destroy_plan(plan%planptr)
        else
          call fftw_destroy_plan(plan%planptr)
        end if
#else
        call fftw_destroy_plan(plan%planptr)
#endif
      end if

      plan%planptr = c_null_ptr
    end subroutine PoisFFT_Plan2D_Finalize

    subroutine PoisFFT_Plan3D_Finalize(plan)
      type(PoisFFT_Plan3D) :: plan
      
      if (c_associated(plan%planptr).and.plan%planowner) then
#ifdef MPI
        if (plan%distributed) then
          call pfft_destroy_plan(plan%planptr)
        else
          call fftw_destroy_plan(plan%planptr)
        end if
#else
        call fftw_destroy_plan(plan%planptr)
#endif
      end if

      plan%planptr = c_null_ptr
    end subroutine PoisFFT_Plan3D_Finalize


    subroutine PoisFFT_InitThreads(nthreads)  !instructs fftw to plan to use nthreads threads
      integer, intent(in) :: nthreads
#if defined(_OPENMP) || defined(ENABLE_PTHREADS)
      integer(c_int) :: error
      error =  fftw_init_threads()

      if (error==0) then
        write(*,*) "Error when initializing FFTW for threads."
      else
        call fftw_plan_with_nthreads(int(nthreads,c_int))
      end if
#endif
    end subroutine PoisFFT_InitThreads

#ifdef MPI
    subroutine PoisFFT_PFFT_Init
      call pfft_init()
    end subroutine

    subroutine PoisFFT_PFFT_InitThreads(nthreads)  !instructs PFFT and fftw to plan to use nthreads threads
      integer, intent(in) :: nthreads
#if defined(_OPENMP) || defined(ENABLE_PTHREADS)
        call pfft_plan_with_nthreads(int(nthreads,c_int))
#endif
    end subroutine PoisFFT_PFFT_InitThreads
#endif

#undef RP
#undef CP
#undef fftw_execute_complex
#undef fftw_execute_real
#undef pfft_execute_gen
#undef fftw_mpi_execute_complex