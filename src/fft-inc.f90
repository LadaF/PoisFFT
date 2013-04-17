  use iso_c_binding
  use PoisFFT_Precisions

  implicit none

  include 'customfftw3.f90'


  integer, parameter :: FFT_Complex = 0
  integer, parameter :: FFT_RealEven00  = FFTW_REDFT00
  integer, parameter :: FFT_RealEven01  = FFTW_REDFT01
  integer, parameter :: FFT_RealEven10  = FFTW_REDFT10
  integer, parameter :: FFT_RealEven11  = FFTW_REDFT11
  integer, parameter :: FFT_RealOdd00   = FFTW_RODFT00
  integer, parameter :: FFT_RealOdd01   = FFTW_RODFT01
  integer, parameter :: FFT_RealOdd10   = FFTW_RODFT10
  integer, parameter :: FFT_RealOdd11   = FFTW_RODFT11

  logical, parameter :: BLOCK_DECOMP = .false.

  type PoisFFT_Plan1D
    type(c_ptr)         :: planptr=c_null_ptr
    logical             :: planowner=.false.
    integer             :: plantype
    integer(c_int)      :: dir
    integer             :: nx
    integer(c_long)     :: cnt
  end type PoisFFT_Plan1D


  type PoisFFT_Plan2D
    type(c_ptr)         :: planptr=c_null_ptr
    logical             :: planowner=.false.
    integer             :: plantype
    integer(c_int)      :: dir
    integer             :: nx
    integer             :: ny
    integer(c_long)     :: cnt
  end type PoisFFT_Plan2D


  type PoisFFT_Plan3D
    type(c_ptr)         :: planptr=c_null_ptr
    logical             :: planowner=.false.
    integer             :: plantype
    integer(c_int)      :: dir
    integer             :: nx
    integer             :: ny
    integer             :: nz
    integer(c_long)     :: cnt
  end type PoisFFT_Plan3D


  type Slice1D
    integer is, ie
  end type Slice1D

  type Slice2D
    integer is, ie, js, je
  end type Slice2D

  type Slice3D
    integer is, ie, js, je, ks, ke
  end type Slice3D

  type PoisFFT_Solver1D
    real(RP) :: dx
    integer  :: nx
    integer  :: cnt
    integer, dimension(2) :: BCs
    type(PoisFFT_Plan1D) :: forward, backward
    complex(CP), dimension(:),&!contiguous,
                       pointer :: cwork => null()
    real(RP), dimension(:),&!contiguous,
                       pointer :: rwork => null()
  end type PoisFFT_Solver1D

  type PoisFFT_Solver2D
    real(RP) :: dx, dy
    integer  :: nx, ny
    integer  :: cnt
    integer, dimension(4) :: BCs
    type(PoisFFT_Plan2D) :: forward, backward
    complex(CP), dimension(:,:),&!contiguous,
                       pointer :: cwork => null()
    real(RP), dimension(:,:),&!contiguous,
                       pointer :: rwork => null()
  end type PoisFFT_Solver2D

  type PoisFFT_Solver3D
    real(RP) :: dx, dy, dz
    integer  :: nx, ny, nz
    integer  :: cnt
    integer, dimension(6) :: BCs
    type(PoisFFT_Plan3D) :: forward, backward
    complex(CP), dimension(:,:,:),&!contiguous,
                       pointer :: cwork => null()
    real(RP), dimension(:,:,:),&!contiguous,
                       pointer :: rwork => null()
    integer :: nthreads = 1
    !will be used in splitting for some boundary conditions
    type(PoisFFT_Solver1D),dimension(:),allocatable :: Solvers1D
    type(PoisFFT_Solver2D),dimension(:),allocatable :: Solvers2D
  end type PoisFFT_Solver3D

  interface assignment(=)
    module procedure PoisFFT_Solver1D_Assign
    module procedure PoisFFT_Solver2D_Assign
    module procedure PoisFFT_Solver3D_Assign
  end interface

  interface assignment(=)
    module procedure PoisFFT_Plan1D_Assign
    module procedure PoisFFT_Plan2D_Assign
    module procedure PoisFFT_Plan3D_Assign
  end interface

  interface data_deallocate
    module procedure data_deallocate_1D_complex
    module procedure data_deallocate_1D_real
    module procedure data_deallocate_2D_complex
    module procedure data_deallocate_2D_real
    module procedure data_deallocate_3D_complex
    module procedure data_deallocate_3D_real
  end interface

  interface data_allocate_real
    module procedure data_allocate_1D_Real
    module procedure data_allocate_2D_Real
    module procedure data_allocate_3D_Real
  end interface

  interface data_allocate_complex
    module procedure data_allocate_1D_Complex
    module procedure data_allocate_2D_Complex
    module procedure data_allocate_3D_Complex
  end interface

  interface Execute
    module procedure PoisFFT_Plan1D_Execute_Real_s
    module procedure PoisFFT_Plan1D_Execute_Complex_s
    module procedure PoisFFT_Plan2D_Execute_Real_s
    module procedure PoisFFT_Plan2D_Execute_Complex_s
    module procedure PoisFFT_Plan3D_Execute_Real_s
    module procedure PoisFFT_Plan3D_Execute_Complex_s
    module procedure PoisFFT_Plan1D_Execute_Real_d
    module procedure PoisFFT_Plan1D_Execute_Complex_d
    module procedure PoisFFT_Plan2D_Execute_Real_d
    module procedure PoisFFT_Plan2D_Execute_Complex_d
    module procedure PoisFFT_Plan3D_Execute_Real_d
    module procedure PoisFFT_Plan3D_Execute_Complex_d
  end interface


  interface Destroy
    module procedure PoisFFT_Plan1D_Destroy
    module procedure PoisFFT_Plan2D_Destroy
    module procedure PoisFFT_Plan3D_Destroy
  end interface




#define cprec CP
#define rprec RP


  contains

    !1D FFT

    subroutine PoisFFT_Solver1D_Assign(DL,DR)
      type(PoisFFT_Solver1D), intent(out) :: DL
      type(PoisFFT_Solver1D), intent(in)  :: DR

      DL%dx = DR%dx

      DL%nx = DR%nx

      DL%cnt = DR%cnt

      DL%BCs = DR%BCs

      DL%forward = DR%forward
      DL%backward = DR%backward
      DL%rwork => DR%rwork
      DL%cwork => DR%cwork
    end subroutine PoisFFT_Solver1D_Assign

    subroutine PoisFFT_Solver2D_Assign(DL,DR)
      type(PoisFFT_Solver2D), intent(out) :: DL
      type(PoisFFT_Solver2D), intent(in)  :: DR

      DL%dx = DR%dx
      DL%dy = DR%dy

      DL%nx = DR%nx
      DL%ny = DR%ny

      DL%cnt = DR%cnt

      DL%BCs = DR%BCs

      DL%forward = DR%forward
      DL%backward = DR%backward
      DL%rwork => DR%rwork
      DL%cwork => DR%cwork
    end subroutine PoisFFT_Solver2D_Assign

    subroutine PoisFFT_Solver3D_Assign(DL,DR)
      type(PoisFFT_Solver3D), intent(out) :: DL
      type(PoisFFT_Solver3D), intent(in)  :: DR

      DL%dx = DR%dx
      DL%dy = DR%dy
      DL%dz = DR%dz

      DL%nx = DR%nx
      DL%ny = DR%ny
      DL%nz = DR%nz

      DL%cnt = DR%cnt

      DL%BCs = DR%BCs

      DL%forward = DR%forward
      DL%backward = DR%backward
      DL%rwork => DR%rwork
      DL%cwork => DR%cwork

      if (allocated(DR%Solvers1D)) then
         allocate(DL%Solvers1D(size(DR%Solvers1D)))
         DL%Solvers1D = DR%Solvers1D
      endif

      if (allocated(DR%Solvers2D)) then
         allocate(DL%Solvers2D(size(DR%Solvers2D)))
         DL%Solvers2D = DR%Solvers2D
      endif

    end subroutine PoisFFT_Solver3D_Assign



    subroutine PoisFFT_Plan1D_Assign(PL,PR)
      type(PoisFFT_Plan1D), intent(out) :: PL
      type(PoisFFT_Plan1D), intent(in)  :: PR

      PL%planptr = PR%planptr
      PL%dir = PL%dir
      PL%nx = PL%nx
      PL%cnt = PL%cnt
    end subroutine PoisFFT_Plan1D_Assign

    subroutine PoisFFT_Plan2D_Assign(PL,PR)
      type(PoisFFT_Plan2D), intent(out) :: PL
      type(PoisFFT_Plan2D), intent(in)  :: PR

      PL%planptr = PR%planptr
      PL%dir = PL%dir
      PL%nx = PL%nx
      PL%ny = PL%ny
      PL%cnt = PL%cnt
    end subroutine PoisFFT_Plan2D_Assign

    subroutine PoisFFT_Plan3D_Assign(PL,PR)
      type(PoisFFT_Plan3d), intent(out) :: PL
      type(PoisFFT_Plan3d), intent(in)  :: PR

      PL%planptr = PR%planptr
      PL%dir = PL%dir
      PL%nx = PL%nx
      PL%ny = PL%ny
      PL%nz = PL%nz
      PL%cnt = PL%cnt
    end subroutine PoisFFT_Plan3D_Assign





    function PoisFFT_Plan1D_QuickCreate(D, plantypes) result(plan)
      type(PoisFFT_Plan1D) :: plan
      type(PoisFFT_Solver1D) :: D
      integer(c_int),intent(in) :: plantypes(:)

      plan = PoisFFT_Plan1D_Create(D, plantypes,&
             Slice1D(1,D%nx),&
             .false.)
    end function PoisFFT_Plan1D_QuickCreate

    function PoisFFT_Plan2D_QuickCreate(D, plantypes) result(plan)
      type(PoisFFT_Plan2D) :: plan
      type(PoisFFT_Solver2D) :: D
      integer(c_int),intent(in) :: plantypes(:)

      plan = PoisFFT_Plan2D_Create(D, plantypes,&
             Slice2D(1,D%nx,1,D%ny),&
             .false.)
    end function PoisFFT_Plan2D_QuickCreate

    function PoisFFT_Plan3D_QuickCreate(D, plantypes) result(plan)
      type(PoisFFT_Plan3D) :: plan
      type(PoisFFT_Solver3D) :: D
      integer(c_int),intent(in) :: plantypes(:)

      plan = PoisFFT_Plan3D_Create(D, plantypes,&
             Slice3D(1,D%nx,1,D%ny,1,D%nz),&
             .false.)
    end function PoisFFT_Plan3D_QuickCreate





    function PoisFFT_Plan1D_Create(D, plantypes, sl, usework) result(plan)
#define dimensions 1
#include "plan_create.f90"
#undef dimensions
    end function PoisFFT_Plan1D_Create

    function PoisFFT_Plan2D_Create(D, plantypes, sl, usework) result(plan)
#define dimensions 2
#include "plan_create.f90"
#undef dimensions
    end function PoisFFT_Plan2D_Create

    function PoisFFT_Plan3D_Create(D, plantypes, sl, usework) result(plan)
#define dimensions 3
#include "plan_create.f90"
#undef dimensions
    end function PoisFFT_Plan3D_Create




    subroutine data_allocate_1D_Complex(D, data)
#define dimensions 1
#define realcomplex 2

#include "data_allocate.f90"

#undef dimensions
#undef realcomplex
    end subroutine data_allocate_1D_Complex

    subroutine data_allocate_1D_Real(D, data)
#define dimensions 1
#define realcomplex 1

#include "data_allocate.f90"

#undef dimensions
#undef realcomplex
    end subroutine data_allocate_1D_Real


    subroutine data_allocate_2D_Complex(D, data)
#define dimensions 2
#define realcomplex 2

#include "data_allocate.f90"

#undef dimensions
#undef realcomplex
    end subroutine data_allocate_2D_Complex

    subroutine data_allocate_2D_Real(D, data)
#define dimensions 2
#define realcomplex 1

#include "data_allocate.f90"

#undef dimensions
#undef realcomplex
    end subroutine data_allocate_2D_Real

    subroutine data_allocate_3D_Complex(D, data)
#define dimensions 3
#define realcomplex 2

#include "data_allocate.f90"

#undef dimensions
#undef realcomplex
    end subroutine data_allocate_3D_Complex

    subroutine data_allocate_3D_Real(D, data)
#define dimensions 3
#define realcomplex 1

#include "data_allocate.f90"

#undef dimensions
#undef realcomplex
    end subroutine data_allocate_3D_Real




    subroutine PoisFFT_Plan1D_Execute_Complex_s(plan, data)
      type(PoisFFT_Plan1D)                   :: plan
      complex(SCP), dimension(:),&!contiguous,
                                     pointer :: data

      call fftwf_execute_dft(plan%planptr, data, data)

    end subroutine PoisFFT_Plan1D_Execute_Complex_s


    subroutine PoisFFT_Plan1D_Execute_Real_s(plan, data)
      type(PoisFFT_Plan1D)                   :: plan
      real(SRP), dimension(:),&!contiguous,
                                     pointer :: data

      call fftwf_execute_r2r(plan%planptr, data, data)

    end subroutine PoisFFT_Plan1D_Execute_Real_s

    subroutine PoisFFT_Plan2D_Execute_Complex_s(plan, data)
      type(PoisFFT_Plan2D)                   :: plan
      complex(SCP), dimension(:,:),&!contiguous,
                                     pointer :: data

      call fftwf_execute_dft(plan%planptr, data, data)

    end subroutine PoisFFT_Plan2D_Execute_Complex_s


    subroutine PoisFFT_Plan2D_Execute_Real_s(plan, data)
      type(PoisFFT_Plan2D)                   :: plan
      real(SRP), dimension(:,:),&!contiguous,
                                     pointer :: data

      call fftwf_execute_r2r(plan%planptr, data, data)

    end subroutine PoisFFT_Plan2D_Execute_Real_s

    subroutine PoisFFT_Plan3D_Execute_Complex_s(plan, data)
      type(PoisFFT_Plan3D)                   :: plan
      complex(SCP), dimension(:,:,:),&!contiguous,
                                     pointer :: data

      call fftwf_execute_dft(plan%planptr, data, data)

    end subroutine PoisFFT_Plan3D_Execute_Complex_s


    subroutine PoisFFT_Plan3D_Execute_Real_s(plan, data)
      type(PoisFFT_Plan3D)                   :: plan
      real(SRP), dimension(:,:,:),&!contiguous,
                                     pointer :: data

      call fftwf_execute_r2r(plan%planptr, data, data)

    end subroutine PoisFFT_Plan3D_Execute_Real_s


    subroutine PoisFFT_Plan1D_Execute_Complex_d(plan, data)
      type(PoisFFT_Plan1D)                   :: plan
      complex(DCP), dimension(:),&!contiguous,
                                     pointer :: data

      call fftw_execute_dft(plan%planptr, data, data)

    end subroutine PoisFFT_Plan1D_Execute_Complex_d


    subroutine PoisFFT_Plan1D_Execute_Real_d(plan, data)
      type(PoisFFT_Plan1D)                   :: plan
      real(DRP), dimension(:),&!contiguous,
                                     pointer :: data

      call fftw_execute_r2r(plan%planptr, data, data)

    end subroutine PoisFFT_Plan1D_Execute_Real_d

    subroutine PoisFFT_Plan2D_Execute_Complex_d(plan, data)
      type(PoisFFT_Plan2D)                   :: plan
      complex(DCP), dimension(:,:),&!contiguous,
                                     pointer :: data

      call fftw_execute_dft(plan%planptr, data, data)

    end subroutine PoisFFT_Plan2D_Execute_Complex_d


    subroutine PoisFFT_Plan2D_Execute_Real_d(plan, data)
      type(PoisFFT_Plan2D)                   :: plan
      real(DRP), dimension(:,:),&!contiguous,
                                     pointer :: data

      call fftw_execute_r2r(plan%planptr, data, data)

    end subroutine PoisFFT_Plan2D_Execute_Real_d

    subroutine PoisFFT_Plan3D_Execute_Complex_d(plan, data)
      type(PoisFFT_Plan3D)                   :: plan
      complex(DCP), dimension(:,:,:),&!contiguous,
                                     pointer :: data

      call fftw_execute_dft(plan%planptr, data, data)

    end subroutine PoisFFT_Plan3D_Execute_Complex_d


    subroutine PoisFFT_Plan3D_Execute_Real_d(plan, data)
      type(PoisFFT_Plan3D)                   :: plan
      real(DRP), dimension(:,:,:),&!contiguous,
                                     pointer :: data

      call fftw_execute_r2r(plan%planptr, data, data)

    end subroutine PoisFFT_Plan3D_Execute_Real_d






    subroutine data_deallocate_1D_complex(data)
      complex(CP), dimension(:), pointer :: data
      type(c_ptr) :: p

      p = c_loc(data(1))
      call fftw_free(p)
      nullify(data)
    end subroutine data_deallocate_1D_complex

    subroutine data_deallocate_1D_real(data)
      real(RP), dimension(:), pointer :: data
      type(c_ptr) :: p

      p = c_loc(data(1))
      call fftw_free(p)
      nullify(data)
    end subroutine data_deallocate_1D_real

    subroutine data_deallocate_2D_complex(data)
      complex(CP), dimension(:,:), pointer :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine data_deallocate_2D_complex

    subroutine data_deallocate_2D_real(data)
      real(RP), dimension(:,:), pointer :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine data_deallocate_2D_real

    subroutine data_deallocate_3D_complex(data)
      complex(CP), dimension(:,:,:), pointer :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine data_deallocate_3D_complex


    subroutine data_deallocate_3D_real(data)
      real(RP), dimension(:,:,:), pointer :: data
      type(c_ptr) :: p

      p = c_loc(data(1,1,1))
      call fftw_free(p)
      nullify(data)
    end subroutine data_deallocate_3D_real





    subroutine PoisFFT_Plan1D_Destroy(plan)
      type(PoisFFT_Plan1D) :: plan

      if (c_associated(plan%planptr).and.plan%planowner)  call fftw_destroy_plan(plan%planptr)

    end subroutine PoisFFT_Plan1D_Destroy

    subroutine PoisFFT_Plan2D_Destroy(plan)
      type(PoisFFT_Plan2D) :: plan

      if (c_associated(plan%planptr).and.plan%planowner)  call fftw_destroy_plan(plan%planptr)

    end subroutine PoisFFT_Plan2D_Destroy

    subroutine PoisFFT_Plan3D_Destroy(plan)
      type(PoisFFT_Plan3D) :: plan

      if (c_associated(plan%planptr).and.plan%planowner)  call fftw_destroy_plan(plan%planptr)

    end subroutine PoisFFT_Plan3D_Destroy


    subroutine PoisFFT_InitThreads(nthreads)  !instructs fftw to plan to use nthreads threads
      integer, intent(in) :: nthreads
      integer(c_int) :: error

      error =  fftw_init_threads()

      if (error==0) then
        write(*,*) "Error when initializing FFTW for threads."
      else
        call fftw_plan_with_nthreads(int(nthreads,c_int))
      end if

    end subroutine PoisFFT_InitThreads

    