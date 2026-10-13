#if (PREC==2)

#define RP DRP
#define CP DCP
#define MPI_RP MPI_DOUBLE_PRECISION
#define MPI_CP MPI_DOUBLE_COMPLEX

#define _RP _DRP

#else

#define RP SRP
#define CP SCP
#define MPI_RP MPI_REAL
#define MPI_CP MPI_COMPLEX

#define _RP _SRP
#define _CP _SCP

#endif

  use iso_c_binding
  !$ use omp_lib
  use PoisFFT_Precisions
  use PoisFFT_Parameters
#ifdef MPI
  use mpi
#endif
  implicit none

  private
  public :: PoisFFT_Solver1D, PoisFFT_Solver2D, PoisFFT_Solver3D, &
            PoisFFT_Solver3D_nonuniform_z, &
            Finalize, Execute

  real(RP), parameter, private :: pi = 3.141592653589793238462_RP

  !result of real_transform_type_forward/backward for an unsupported pair
  !of boundary conditions, differs from all FFTW kinds and FFT_Complex
  integer, parameter, private :: FFT_Unsupported = -2


  interface PoisFFT_Solver1D
    module procedure PoisFFT_Solver1D__New
  end interface

  interface PoisFFT_Solver2D
    module procedure PoisFFT_Solver2D__New
  end interface

  interface PoisFFT_Solver3D
    module procedure PoisFFT_Solver3D__New
  end interface

  interface PoisFFT_Solver3D_nonuniform_z
    module procedure PoisFFT_Solver3D_nonuniform_z__New
  end interface

  interface Finalize
    module procedure PoisFFT_Solver1D__Finalize
    module procedure PoisFFT_Solver2D__Finalize
    module procedure PoisFFT_Solver3D__Finalize
    module procedure PoisFFT_Solver3D_nonuniform_z__Finalize
  end interface Finalize


  interface Init
    module procedure PoisFFT_Solver1D_Init
    module procedure PoisFFT_Solver2D_Init
    module procedure PoisFFT_Solver3D_Init
    module procedure PoisFFT_Solver3D_nonuniform_z_Init
  end interface Init

  interface Execute
    module procedure PoisFFT_Solver1D__Execute
    module procedure PoisFFT_Solver2D__Execute
    module procedure PoisFFT_Solver3D__Execute
    module procedure PoisFFT_Solver3D_nonuniform_z__Execute
  end interface Execute


  contains


    function PoisFFT_Solver3D__New(nxyz,Lxyz,BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads,planning) result(D)
      type(PoisFFT_Solver3D) :: D

      integer, intent(in)   :: nxyz(3)
      real(RP), intent(in)  :: Lxyz(3)
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(3)
      integer, intent(in), optional :: offs(3)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads
      integer, intent(in), optional :: planning

      D%nxyz = nxyz

      D%Lx = Lxyz(1)
      D%Ly = Lxyz(2)
      D%Lz = Lxyz(3)

      D%nx = nxyz(1)
      D%ny = nxyz(2)
      D%nz = nxyz(3)
      
      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
        D%gny = gnxyz(2)
        D%gnz = gnxyz(3)
      else
        D%gnx = D%nx
        D%gny = D%ny
        D%gnz = D%nz
      end if

      if (present(offs)) then
        D%offx = offs(1)
        D%offy = offs(2)
        D%offz = offs(3)
      end if

      D%cnt = product(D%nxyz)
      
      D%gcnt = product(int([D%gnx, D%gny, D%gnz], kind(D%gcnt)))

      D%BCs = BCs
      
      if (present(approximation)) D%approximation = approximation
      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
#ifdef MPI      
      else
        stop "No PFFT comm present in PoisFFT_Solver3D__New."
#endif
      end if

      if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      if (present(planning)) D%planning = planning

      !create fftw plans and allocate working arrays
      call Init(D)
    end function PoisFFT_Solver3D__New
    
    


    subroutine PoisFFT_Solver3D_Init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: real_forw, real_back
      integer :: real_forw_xy, real_back_xy
      integer :: real_forw_yz, real_back_yz
      integer :: real_forw_z, real_back_z
      integer :: real_forw_x, real_back_x
      integer :: i
      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel

      call check_BC_pairs(D%BCs)

      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2)) * &
                      norm_factor(D%gny,D%BCs(3:4)) * &
                      norm_factor(D%gnz,D%BCs(5:6))
      
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))
      allocate(D%denomz(D%nz))
      
      
#ifdef MPI
      call PoisFFT_PFFT_init()
#else
      if (D%nthreads>1) call PoisFFT_InitThreads(1)
#endif

      if (all(D%BCs==PoisFFT_Periodic)) then

#ifdef MPI
        if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)
#else
        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)
#endif

        call allocate_fftw_complex(D)

        D%forward = PoisFFT_Plan3D(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan3D(D, [FFT_Complex, FFTW_BACKWARD])
        
      else if (all(D%BCs==PoisFFT_Dirichlet) .or. &
               all(D%BCs==PoisFFT_DirichletStag) .or. &
               all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

#ifdef MPI
        if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)
#else
        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)
#endif

        call allocate_fftw_real(D)
        
        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))

        D%forward = PoisFFT_Plan3D(D, [(real_forw, i=1,3)])
        D%backward = PoisFFT_Plan3D(D, [(real_back, i=1,3)])

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag) )) then

#ifdef MPI
        allocate(D%Solvers1D(3:2+D%nthreads))

        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(3))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw, real_forw])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back, real_back])

        do i = 4, 2+D%nthreads
          D%Solvers1D(i) = D%Solvers1D(3)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        if (D%ny<D%gny) then
          if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)

          allocate(D%Solvers2D(1))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

          
        else
          allocate(D%Solvers2D(D%nthreads))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                  [FFT_Complex, FFTW_FORWARD], &
                                                  distributed=.false.)
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                   [FFT_Complex, FFTW_BACKWARD], &
                                                   distributed=.false.)

          do i = 2, D%nthreads
            D%Solvers2D(i) = D%Solvers2D(1)
            D%Solvers2D(i)%forward%planowner = .false.
            D%Solvers2D(i)%backward%planowner = .false.

            call allocate_fftw_complex(D%Solvers2D(i))
          end do
        end if

        call Init_MPI_Buffers(D, 3)

#else

        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(1))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
        
!MPI
#endif

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
                ( (any(D%BCs(5:6)==PoisFFT_DirichletStag) .or.any(D%BCs(5:6)==PoisFFT_Dirichlet)) .and. &
                  all(D%BCs(5:6)==PoisFFT_DirichletStag .or. &
                      D%BCs(5:6)==PoisFFT_Dirichlet .or. &
                      D%BCs(5:6)==PoisFFT_NeumannStag .or. &
                      D%BCs(5:6)==PoisFFT_Neumann))) then

#ifdef MPI
        allocate(D%Solvers1D(3:2+D%nthreads))

        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(3))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw, real_forw])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back, real_back])

        do i = 4, 2+D%nthreads
          D%Solvers1D(i) = D%Solvers1D(3)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        if (D%ny<D%gny) then
          if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)

          allocate(D%Solvers2D(1))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

          
        else
          allocate(D%Solvers2D(D%nthreads))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                  [FFT_Complex, FFTW_FORWARD], &
                                                  distributed=.false.)
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                   [FFT_Complex, FFTW_BACKWARD], &
                                                   distributed=.false.)

          do i = 2, D%nthreads
            D%Solvers2D(i) = D%Solvers2D(1)
            D%Solvers2D(i)%forward%planowner = .false.
            D%Solvers2D(i)%backward%planowner = .false.

            call allocate_fftw_complex(D%Solvers2D(i))
          end do
        end if

        call Init_MPI_Buffers(D, 3)

#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(1))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
      
       !MPI
#endif


      else if (all(D%BCs(3:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(1:2)==PoisFFT_Neumann) .or. &
                all(D%BCs(1:2)==PoisFFT_NeumannStag))) then

#ifdef MPI
        call PoisFFT_Solver3D_Separable_MPI_Init(D)
#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)

        call allocate_fftw_real(D%Solvers1D(1))

        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,1)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
#endif

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:4)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:4)==PoisFFT_NeumannStag) )) then

#ifdef MPI
        call PoisFFT_Solver3D_Separable_MPI_Init(D)
#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,2)

        call allocate_fftw_real(D%Solvers1D(1))

        real_forw = real_transform_type_forward(D%BCs(3:4))
        real_back = real_transform_type_backward(D%BCs(3:4))
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw, real_forw])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back, real_back])

        do i = 2, D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,2)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do
#endif

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:6)==PoisFFT_NeumannStag) )) then

        real_forw = real_transform_type_forward(D%BCs(3:4))
        real_back = real_transform_type_backward(D%BCs(3:4))
        
#ifdef MPI
        !1..x, 2..y, 3..z, not different threads, but directions
        allocate(D%Solvers1D(3*D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)
        D%Solvers1D(2) = PoisFFT_Solver1D_From3D(D,2)
        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)
        
        call allocate_fftw_complex(D%Solvers1D(1))
        call allocate_fftw_real(D%Solvers1D(2))
        call allocate_fftw_real(D%Solvers1D(3))
        
        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_BACKWARD])
        D%Solvers1D(2)%forward = PoisFFT_Plan1D(D%Solvers1D(2), [real_forw])
        D%Solvers1D(2)%backward = PoisFFT_Plan1D(D%Solvers1D(2), [real_back])
        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back])
        
        do i = 4, 1 + 3*(D%nthreads-1), 3
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_complex(D%Solvers1D(i))
        end do

        do i = 5, 2 + 3*(D%nthreads-1), 3
          D%Solvers1D(i) = D%Solvers1D(2)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do
        do i = 6, 3 + 3*(D%nthreads-1), 3
          D%Solvers1D(i) = D%Solvers1D(3)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        call Init_MPI_Buffers(D, 2)
        call Init_MPI_Buffers(D, 3)

#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)

        call allocate_fftw_complex(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_BACKWARD])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_complex(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,1)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw, real_forw])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back, real_back])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
#endif




      else if ((all(D%BCs(1:4)==PoisFFT_Dirichlet) .or. &
                all(D%BCs(1:4)==PoisFFT_DirichletStag)) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic)) then

        real_forw = real_transform_type_forward(D%BCs(1:2))
        real_back = real_transform_type_backward(D%BCs(1:2))
        
#ifdef MPI
        call PoisFFT_Solver3D_Separable_MPI_Init(D)
#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_Complex, FFTW_BACKWARD])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_complex(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw, real_forw])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back, real_back])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
#endif

      else if ( (all(D%BCs==PoisFFT_Neumann .or. &
                     D%BCs==PoisFFT_NeumannStag .or. &
                     D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
               .and. & 
               (all(D%BCs(1:2)==D%BCs(3:4)) .and. &
                any(D%BCs(1:2)/=D%BCs(5:6))) &
               .and. &
               (any(D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
              ) then

        real_forw_xy = real_transform_type_forward(D%BCs(1:2))
        real_back_xy = real_transform_type_backward(D%BCs(1:2))
        real_forw_z = real_transform_type_forward(D%BCs(5:6))
        real_back_z = real_transform_type_backward(D%BCs(5:6))
        
#ifdef MPI
        call PoisFFT_Solver3D_Separable_MPI_Init(D)
#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw_z, real_forw_z])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back_z, real_back_z])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw_xy, real_forw_xy])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back_xy, real_back_xy])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
#endif

      else if ( (all(D%BCs==PoisFFT_Neumann .or. &
                     D%BCs==PoisFFT_NeumannStag .or. &
                     D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
               .and. & 
               (all(D%BCs(3:4)==D%BCs(5:6)) .and. &
                any(D%BCs(1:2)/=D%BCs(3:4))) &
               .and. &
               (any(D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
              ) then

        real_forw_yz = real_transform_type_forward(D%BCs(3:4))
        real_back_yz = real_transform_type_backward(D%BCs(3:4))
        real_forw_x = real_transform_type_forward(D%BCs(1:2))
        real_back_x = real_transform_type_backward(D%BCs(1:2))
        
#ifdef MPI
        call PoisFFT_Solver3D_Separable_MPI_Init(D)
#else
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,1)

        call allocate_fftw_real(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [real_forw_x, real_forw_x])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [real_back_x, real_back_x])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,1)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [real_forw_yz, real_forw_yz])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [real_back_yz, real_back_yz])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
#endif


      else
        !any other combination, every pair was checked by check_BC_pairs

#ifdef MPI
        call PoisFFT_Solver3D_Separable_MPI_Init(D)
#else
        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call allocate_fftw_real(D)

        !the kinds in the C order of dimensions: z, y, x
        D%forward = PoisFFT_Plan3D(D, [real_transform_type_forward(D%BCs(5:6)), &
                                       real_transform_type_forward(D%BCs(3:4)), &
                                       real_transform_type_forward(D%BCs(1:2))])
        D%backward = PoisFFT_Plan3D(D, [real_transform_type_backward(D%BCs(5:6)), &
                                        real_transform_type_backward(D%BCs(3:4)), &
                                        real_transform_type_backward(D%BCs(1:2))])
#endif

      endif


      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD2, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues(eig_fn_FD2, D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD4, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues(eig_fn_FD4, D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_spectral, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues(eig_fn_spectral, D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      end if
      
    end subroutine PoisFFT_Solver3D_Init
    
    


    subroutine PoisFFT_Solver3D__Execute(D,Phi,RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)

      integer   :: ngPhi(3), ngRHS(3)

      ngPhi = (ubound(Phi)-[D%nx,D%ny,D%nz])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny,D%nz])/2

      if (all(D%BCs==PoisFFT_Periodic)) then

        call PoisFFT_Solver3D_FullPeriodic(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs==PoisFFT_Dirichlet) .or. &
               all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver3D_FullDirichlet(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver3D_FullNeumann(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag))) then
        call PoisFFT_Solver3D_PPNs(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))


      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
                ( (any(D%BCs(5:6)==PoisFFT_DirichletStag) .or.any(D%BCs(5:6)==PoisFFT_Dirichlet)) .and. &
                  all(D%BCs(5:6)==PoisFFT_DirichletStag .or. &
                      D%BCs(5:6)==PoisFFT_Dirichlet .or. &
                      D%BCs(5:6)==PoisFFT_NeumannStag .or. &
                      D%BCs(5:6)==PoisFFT_Neumann))) then
        call PoisFFT_Solver3D_PPD(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(3:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(1:2)==PoisFFT_Neumann) .or. &
                all(D%BCs(1:2)==PoisFFT_NeumannStag))) then
        call PoisFFT_Solver3D_NsPP(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:6)==PoisFFT_NeumannStag))) then

        call PoisFFT_Solver3D_PNsNs(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:4)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:4)==PoisFFT_NeumannStag) )) then
                
        call PoisFFT_Solver3D_PNsP(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if ((all(D%BCs(1:4)==PoisFFT_Dirichlet) .or. &
                all(D%BCs(1:4)==PoisFFT_DirichletStag)) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic)) then

        call PoisFFT_Solver3D_DDP(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))
                     
      else if ( (all(D%BCs==PoisFFT_Neumann .or. &
                     D%BCs==PoisFFT_NeumannStag .or. &
                     D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
               .and. & 
               (all(D%BCs(1:2)==D%BCs(3:4)) .and. &
                any(D%BCs(1:2)/=D%BCs(5:6))) &
               .and. &
               (any(D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
              ) then

        call PoisFFT_Solver3D_2real1real(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if ( (all(D%BCs==PoisFFT_Neumann .or. &
                     D%BCs==PoisFFT_NeumannStag .or. &
                     D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
               .and. & 
               (all(D%BCs(3:4)==D%BCs(5:6)) .and. &
                any(D%BCs(1:2)/=D%BCs(3:4))) &
               .and. &
               (any(D%BCs==PoisFFT_Dirichlet .or. &
                     D%BCs==PoisFFT_DirichletStag)) &
              ) then

        call PoisFFT_Solver3D_1real2real(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else

        call PoisFFT_Solver3D_Separable(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      endif

    end subroutine PoisFFT_Solver3D__Execute
    
    
    
    
    subroutine PoisFFT_Solver3D__Finalize(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: i

      call Finalize(D%forward)
      call Finalize(D%backward)

      if (associated(D%rwork)) call deallocate_fftw(D%rwork)
      if (associated(D%cwork)) call deallocate_fftw(D%cwork)

      if (allocated(D%Solvers1D)) then
        do i = lbound(D%Solvers1D,1), ubound(D%Solvers1D,1)
          call Finalize(D%Solvers1D(i))
        end do

        deallocate(D%Solvers1D)
      endif

      if (allocated(D%Solvers2D)) then
        do i = lbound(D%Solvers2D,1), ubound(D%Solvers2D,1)
          call Finalize(D%Solvers2D(i))
        end do

        deallocate(D%Solvers2D)
      endif

    endsubroutine PoisFFT_Solver3D__Finalize







    function PoisFFT_Solver1D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver1D)             :: D
      class(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction
#ifdef MPI
      integer :: ie, dims
      
      call MPI_Cartdim_get(D3D%mpi%comm, dims, ie)
      if (ie/=0) stop "Error executing MPI_Cartdim_get."

      if (dims==2) then
        !We see the dimensions reversed in Fortran!
        if (direction==1) then
          D%mpi%comm = MPI_COMM_SELF
        else if (direction==2) then
          call MPI_Cart_sub(D3D%mpi%comm, [.false.,.true.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
          call MPI_Comm_size(D%mpi%comm, D%mpi%np, ie)
          D%mpi_transpose_needed = D%mpi%np > 1
        else
          call MPI_Cart_sub(D3D%mpi%comm, [.true.,.false.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
          call MPI_Comm_size(D%mpi%comm, D%mpi%np, ie)
          D%mpi_transpose_needed = D%mpi%np > 1
        end if
      else
        stop "Not implemented."
      end if

      call MPI_Comm_rank(D%mpi%comm, D%mpi%rank, ie)

#endif

      if (direction==1) then
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%BCs = D3D%BCs(1:2)
      else if (direction==2) then
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%BCs = D3D%BCs(3:4)
      else
        D%nx = D3D%nz
        D%gnx = D3D%gnz
        D%offx = D3D%offz
        D%BCs = D3D%BCs(5:6)
      endif
      
      D%nxyz = [D%nx]

      D%cnt = D%nx
      
      D%gcnt = int(D%gnx, kind(D%gcnt))

      D%planning = D3D%planning
    end function PoisFFT_Solver1D_From3D



    function PoisFFT_Solver2D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver2D)             :: D
      class(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction
#ifdef MPI
      integer :: ie, dims
      integer :: nps(2), coords(2)
      logical :: pers(2)
      
      call MPI_Cartdim_get(D3D%mpi%comm, dims, ie)
      if (ie/=0) stop "Error executing MPI_Cartdim_get."
        
      if (dims==2) then
        !We see the dimensions reversed in Fortran!
        if (direction==1) then
          call MPI_Cart_get(D3D%mpi%comm, dims, nps, pers, coords, ie)
          if (ie/=0) stop "Error executing MPI_Cart_get."
      
          if (product(nps)==1) then
            D%mpi%np = 1
            D%mpi%comm = MPI_COMM_SELF
            D%mpi%comm_dim = 0
          else if (nps(1)==1) then
            call MPI_Cart_sub(D3D%mpi%comm, [.false.,.true.], D%mpi%comm, ie)
            if (ie/=0) stop "Error executing MPI_Cart_sub."
            D%mpi%comm_dim = 1
          else if (nps(2)==1) then
            call MPI_Cart_sub(D3D%mpi%comm, [.true.,.false.], D%mpi%comm, ie)
            if (ie/=0) stop "Error executing MPI_Cart_sub."
            D%mpi%comm_dim = 1
          else
            D%mpi%comm = D3D%mpi%comm
            D%mpi%comm_dim = 2
          end if
        else if (direction==2) then
          call MPI_Cart_sub(D3D%mpi%comm, [.true.,.false.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
          D%mpi%comm_dim = 1
        else
          call MPI_Cart_sub(D3D%mpi%comm, [.false.,.true.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
          D%mpi%comm_dim = 1
        end if
      else
        stop "Not implemented."
      end if
      
#endif

      if (direction==1) then
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%offy = D3D%offz
        D%BCs = D3D%BCs(3:6)
      else if (direction==2) then
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%offy = D3D%offz
        D%BCs = D3D%BCs([1,2,5,6])
      else
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%ny = D3D%ny
        D%gny = D3D%gny
        D%offy = D3D%offy
        D%BCs = D3D%BCs(1:4)
      endif

      D%nxyz = [D%nx, D%ny]

      D%cnt = D%nx * D%ny
      
      D%gcnt = int(D%gnx, kind(D%gcnt)) * int(D%gny, kind(D%gcnt))

      D%planning = D3D%planning
    end function PoisFFT_Solver2D_From3D




    function PoisFFT_Solver2D__New(nxyz,Lxyz,BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads,planning) result(D)
      type(PoisFFT_Solver2D) :: D

      integer, intent(in)   :: nxyz(2)
      real(RP), intent(in)  :: Lxyz(2)
      integer, intent(in)   :: bcs(4)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(2)
      integer, intent(in), optional :: offs(2)
      integer, intent(in), optional  :: mpi_comm
      integer, intent(in), optional :: nthreads
      integer, intent(in), optional :: planning

      D%nxyz = nxyz

      D%Lx = Lxyz(1)
      D%Ly = Lxyz(2)

      D%nx = nxyz(1)
      D%ny = nxyz(2)

      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
        D%gny = gnxyz(2)
      else
        D%gnx = D%nx
        D%gny = D%ny
      end if

      if (present(offs)) then
        D%offx = offs(1)
        D%offy = offs(2)
      end if

      D%cnt = product(D%nxyz)

      D%gcnt = product(int([D%gnx, D%gny], kind(D%gcnt)))

      D%BCs = BCs

      if (present(approximation)) D%approximation = approximation
      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
#ifdef MPI      
      else
        stop "No PFFT comm present in PoisFFT_Solver2D__New."
#endif
      end if

      if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      if (present(planning)) D%planning = planning

      !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver2D__New




    subroutine Poisfft_Solver2D_Init(D)
      type(PoisFFT_Solver2D), intent(inout) :: D
        
      call check_BC_pairs(D%BCs)

      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2)) * &
                      norm_factor(D%gny,D%BCs(3:4))
      
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call allocate_fftw_complex(D)

        D%forward = PoisFFT_Plan2D(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan2D(D, [FFT_Complex, FFTW_BACKWARD])

      else
        !any other combination, every pair was checked by check_BC_pairs

        call allocate_fftw_real(D)

        !the kinds in the C order of dimensions: y, x
        D%forward = PoisFFT_Plan2D(D, [real_transform_type_forward(D%BCs(3:4)), &
                                       real_transform_type_forward(D%BCs(1:2))])
        D%backward = PoisFFT_Plan2D(D, [real_transform_type_backward(D%BCs(3:4)), &
                                        real_transform_type_backward(D%BCs(1:2))])

      endif
      
      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD2, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD4, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_spectral, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      end if
    end subroutine PoisFFT_Solver2D_Init




    subroutine PoisFFT_Solver2D__Execute(D,Phi,RHS)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:)
      real(RP), intent(in)  :: RHS(:,:)

      integer   :: ngPhi(2), ngRHS(2)


      ngPhi = (ubound(Phi)-[D%nx,D%ny])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny])/2


      if (all(D%BCs==PoisFFT_Periodic)) then

        call PoisFFT_Solver2D_FullPeriodic(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny))

      else

        call PoisFFT_Solver2D_Separable(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny))

      endif

    end subroutine PoisFFT_Solver2D__Execute





    subroutine PoisFFT_Solver2D__Finalize(D)
      type(PoisFFT_Solver2D), intent(inout) :: D

      call Finalize(D%forward)
      call Finalize(D%backward)

      if (associated(D%rwork)) call deallocate_fftw(D%rwork)
      if (associated(D%cwork)) call deallocate_fftw(D%cwork)

    endsubroutine PoisFFT_Solver2D__Finalize






















    function PoisFFT_Solver1D__New(nxyz,Lxyz,BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads,planning) result(D)
      type(PoisFFT_Solver1D) :: D

      integer, intent(in)   :: nxyz(1)
      real(RP), intent(in)  :: Lxyz(1)
      integer, intent(in)   :: bcs(2)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(1)
      integer, intent(in), optional :: offs(1)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads
      integer, intent(in), optional :: planning
#ifdef MPI
      integer :: ie, dims
#endif
      D%nxyz = nxyz
      
      D%Lx = Lxyz(1)

      D%nx = nxyz(1)

      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
      else
        D%gnx = D%nx
      end if

      if (present(offs)) then
        D%offx = offs(1)
      end if

      D%cnt = D%nx

      D%gcnt = int(D%gnx, kind(D%gcnt))

      D%BCs = BCs

      if (present(approximation)) D%approximation = approximation

      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
#ifdef MPI
        call MPI_Cartdim_get(D%mpi%comm, dims, ie)
        if (ie/=0) stop "Error executing MPI_Cartdim_get."
        D%mpi_transpose_needed = dims > 0
      else
        stop "No PFFT comm present in PoisFFT_Solver1D__New."
#endif
      end if
      
      
      
     if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      if (present(planning)) D%planning = planning

       !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver1D__New




    subroutine Poisfft_Solver1D_Init(D)
      type(PoisFFT_Solver1D), intent(inout) :: D
      integer :: real_forw, real_back
      
      call check_BC_pairs(D%BCs)

      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2))
      
      allocate(D%denomx(D%gnx))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call allocate_fftw_complex(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan1D(D, [FFT_Complex, FFTW_BACKWARD])

      else 
      
        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [real_transform_type_forward(D%BCs(1:2))])
        D%backward = PoisFFT_Plan1D(D, [real_transform_type_backward(D%BCs(1:2))])
        
      endif
      
      !for the whole grid, with MPI it is solved on the first process
      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%gnx, D%gnx, 0)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%gnx, D%gnx, 0)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%gnx, D%gnx, 0)
      end if

    end subroutine PoisFFT_Solver1D_Init
    
    
    
    
    subroutine PoisFFT_Solver1D__Execute(D,Phi,RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)

      integer   :: ngPhi(1), ngRHS(1)


      ngPhi = (ubound(Phi)-D%nx)/2
      ngRHS = (ubound(RHS)-D%nx)/2

#ifdef MPI
      if (D%mpi_transpose_needed) then
        call PoisFFT_Solver1D_MPI(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))
        return
      end if
#endif

      call PoisFFT_Solver1D_Local(D,&
               Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
               RHS(ngRHS(1)+1:ngRHS(1)+D%nx))

    end subroutine PoisFFT_Solver1D__Execute


    !Solution of the whole grid in Phi and RHS (without ghost points).
    subroutine PoisFFT_Solver1D_Local(D,Phi,RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)

      if (all(D%BCs==PoisFFT_Periodic)) then

        call PoisFFT_Solver1D_FullPeriodic(D, Phi, RHS)

      else if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver1D_FullNeumann(D, Phi, RHS)

      else
        !Dirichlet on at least one end, every pair was checked by check_BC_pairs

        call PoisFFT_Solver1D_FullDirichlet(D, Phi, RHS)

      endif

    end subroutine PoisFFT_Solver1D_Local


#ifdef MPI
    !The 1D solver with a distributed grid (local parts without ghost
    !points): the right-hand side is gathered on the first process of the
    !communicator, solved there and the solution is scattered back.
    subroutine PoisFFT_Solver1D_MPI(D,Phi,RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      integer, allocatable :: displs(:), counts(:)
      real(RP), allocatable :: grhs(:), gphi(:)
      integer :: ie

      call MPI_Comm_rank(D%mpi%comm, D%mpi%rank, ie)
      call MPI_Comm_size(D%mpi%comm, D%mpi%np, ie)
      allocate(displs(D%mpi%np), counts(D%mpi%np))
      if (D%mpi%rank==0) then
        allocate(grhs(D%gnx), gphi(D%gnx))
      else
        allocate(grhs(1), gphi(1))
      end if

      call MPI_Gather(D%offx, 1, MPI_INTEGER, &
                      displs, 1, MPI_INTEGER, 0, &
                      D%mpi%comm, ie)

      call MPI_Gather(D%nx, 1, MPI_INTEGER, &
                      counts, 1, MPI_INTEGER, 0, &
                      D%mpi%comm, ie)

      call MPI_Gatherv(RHS, D%nx, MPI_RP, &
                       grhs, counts, displs, MPI_RP, 0, &
                       D%mpi%comm, ie)
      if (ie/=0) stop "Error in MPI_Gatherv!"

      if (D%mpi%rank==0) call PoisFFT_Solver1D_Local(D, gphi, grhs)

      call MPI_Scatterv(gphi, counts, displs, MPI_RP, &
                        Phi, D%nx, MPI_RP, 0, &
                        D%mpi%comm, ie)
      if (ie/=0) stop "Error in MPI_Scatterv!"
    end subroutine PoisFFT_Solver1D_MPI
#endif




    subroutine PoisFFT_Solver1D__Finalize(D)
      type(PoisFFT_Solver1D), intent(inout) :: D

      call Finalize(D%forward)
      call Finalize(D%backward)

      if (associated(D%rwork)) call deallocate_fftw(D%rwork)
      if (associated(D%cwork)) call deallocate_fftw(D%cwork)

    endsubroutine PoisFFT_Solver1D__Finalize





    
    
    
    
    
    
    
    
    !If ierr is present, errors in the arguments are returned in ierr
    !(PoisFFT_ErrorBCs, PoisFFT_ErrorGrid, PoisFFT_OK otherwise) and in
    !errmsg, the returned solver can only be passed to Finalize. Without
    !ierr the program stops with the error message.
    function PoisFFT_Solver3D_nonuniform_z__New(nxyz,Lxy,z,z_u, &
                                  BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads,ierr,planning, &
                                  errmsg) result(D)
      type(PoisFFT_Solver3D_nonuniform_z) :: D

      integer, intent(in)   :: nxyz(3)
      real(RP), intent(in)  :: Lxy(2)
      real(RP), intent(in)  :: z(:), z_u(0:)
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(3)
      integer, intent(in), optional :: offs(3)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads
      integer, intent(in), optional :: planning
      integer, intent(out),optional :: ierr
      character(*), intent(out), optional :: errmsg

      if (present(ierr)) ierr = PoisFFT_OK
      if (present(errmsg)) errmsg = ""

      D%nxyz = nxyz

      D%Lx = Lxy(1)
      D%Ly = Lxy(2)
      D%Lz = 0

      D%nx = nxyz(1)
      D%ny = nxyz(2)
      D%nz = nxyz(3)
      
      allocate(D%z(1:ubound(z,1)))
      D%z = z
      allocate(D%z_u(0:ubound(z_u,1)))
      D%z_u(:) = z_u
      
      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
        D%gny = gnxyz(2)
        D%gnz = gnxyz(3)
      else
        D%gnx = D%nx
        D%gny = D%ny
        D%gnz = D%nz
      end if

      if (present(offs)) then
        D%offx = offs(1)
        D%offy = offs(2)
        D%offz = offs(3)
      end if

      D%cnt = product(D%nxyz)
      
      D%gcnt = product(int([D%gnx, D%gny, D%gnz], kind(D%gcnt)))

      D%BCs = BCs
      
      if (present(approximation)) D%approximation = approximation
      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
#ifdef MPI      
      else
        stop "No PFFT comm present in PoisFFT_Solver3D__New."
#endif
      end if

      if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      if (present(planning)) D%planning = planning

      !create fftw plans and allocate working arrays
      call Init(D, ierr, errmsg)
    end function PoisFFT_Solver3D_nonuniform_z__New
    
    
    
    
    subroutine PoisFFT_Solver3D_nonuniform_z_Init(D, ierr, errmsg)
      type(PoisFFT_Solver3D_nonuniform_z), intent(inout) :: D
      integer,      intent(out), optional :: ierr
      character(*), intent(out), optional :: errmsg
      integer :: real_forw, real_back
      integer :: real_forw_xy(2), real_back_xy(2)
      integer :: i, n
      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel

      if (present(ierr)) ierr = PoisFFT_OK
      if (present(errmsg)) errmsg = ""

      !the arguments are checked before anything is allocated, so that
      !a solver returned with an error can be finalized
      !any supported pairs in x and y; z is a staggered grid, the boundaries
      !are at the faces z_u(0) and z_u(n), only staggered conditions in z
      if (.not.BC_pairs_supported(D%BCs) .or. &
          .not.all(D%BCs(5:6)==PoisFFT_NeumannStag .or. &
                   D%BCs(5:6)==PoisFFT_DirichletStag)) then
        call report_error(PoisFFT_ErrorBCs)
        return
      end if

      n = D%gnz
      if (ubound(D%z,1) < n .or. ubound(D%z_u,1) < n) then
        call report_error(PoisFFT_ErrorGrid, 1)
        return
      end if

      if (any(D%z(2:n) <= D%z(1:n-1)) .or. any(D%z_u(1:n) <= D%z_u(0:n-1))) then
        call report_error(PoisFFT_ErrorGrid, 2)
        return
      end if

      !the distance from the point to the face with the Dirichlet condition
      if ((D%BCs(5)==PoisFFT_DirichletStag .and. D%z(1) <= D%z_u(0)) .or. &
          (D%BCs(6)==PoisFFT_DirichletStag .and. D%z(n) >= D%z_u(n))) then
        call report_error(PoisFFT_ErrorGrid, 2)
        return
      end if

      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2)) * &
                      norm_factor(D%gny,D%BCs(3:4))
      
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))
      
      
#ifdef MPI
      call PoisFFT_PFFT_init()
#else
      if (D%nthreads>1) call PoisFFT_InitThreads(1)
#endif

      if (all(D%BCs(1:4)==PoisFFT_Periodic)) then

#ifdef MPI
        call allocate_fftw_complex(D)
        
        allocate(D%Solvers1D(3:2+D%nthreads))

        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(3))

        real_forw = real_transform_type_forward(D%BCs(5:6))
        real_back = real_transform_type_backward(D%BCs(5:6))
        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [real_forw, real_forw])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [real_back, real_back])

        do i = 4, 2+D%nthreads
          D%Solvers1D(i) = D%Solvers1D(3)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do

        if (D%ny<D%gny) then
          if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)

          allocate(D%Solvers2D(1))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D, 3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

          
        else
          allocate(D%Solvers2D(D%nthreads))

          D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

          call allocate_fftw_complex(D%Solvers2D(1))

          D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                  [FFT_Complex, FFTW_FORWARD], &
                                                  distributed=.false.)

          D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), &
                                                   [FFT_Complex, FFTW_BACKWARD], &
                                                   distributed=.false.)

          do i = 2, D%nthreads
            D%Solvers2D(i) = D%Solvers2D(1)
            D%Solvers2D(i)%forward%planowner = .false.
            D%Solvers2D(i)%backward%planowner = .false.

            call allocate_fftw_complex(D%Solvers2D(i))
          end do

        end if

        call Init_MPI_Buffers(D, 3, complex=.true.)

#else

        call allocate_fftw_complex(D)

        allocate(D%Solvers1D(D%nthreads))

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i = 2, D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_complex(D%Solvers2D(i))
        end do

!MPI
#endif

      else
        !any other combination in x and y

#ifdef MPI               
        if (D%nthreads>1) call PoisFFT_PFFT_InitThreads(D%nthreads)
#endif
        !the kinds in the C order of dimensions: y, x
        real_forw_xy = [real_transform_type_forward(D%BCs(3:4)), &
                        real_transform_type_forward(D%BCs(1:2))]
        real_back_xy = [real_transform_type_backward(D%BCs(3:4)), &
                        real_transform_type_backward(D%BCs(1:2))]

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call allocate_fftw_real(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D(D%Solvers2D(1), real_forw_xy)
        D%Solvers2D(1)%backward = PoisFFT_Plan2D(D%Solvers2D(1), real_back_xy)

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.

          call allocate_fftw_real(D%Solvers2D(i))

        end do
        
        allocate(D%Solvers1D(3))
        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

#ifdef MPI
        call Init_MPI_Buffers(D, 3)
#endif

      endif


      if (D%approximation==PoisFFT_FiniteDifference2) then
        D%denomx = eigenvalues(eig_fn_FD2, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD2, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      else if (D%approximation==PoisFFT_FiniteDifference4) then
        D%denomx = eigenvalues(eig_fn_FD4, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_FD4, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      else
        D%denomx = eigenvalues(eig_fn_spectral, D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues(eig_fn_spectral, D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
      end if
      
      allocate(D%mat_a(2:D%gnz))
      allocate(D%mat_b(1:D%gnz))
      allocate(D%mat_c(1:max(1,D%gnz-1)))
    
      do i = 2, D%gnz
        D%mat_a(i) = 1 / (D%z(i)-D%z(i-1)) / (D%z_u(i) - D%z_u(i-1))
      end do
      i = 1
      if (D%gnz==1) then
        D%mat_b(1) = 0
      else
        D%mat_b(1) = -1 / (D%z(i+1)-D%z(i)) / (D%z_u(i) - D%z_u(i-1))
      end if
      do i = 2, D%gnz-1
        D%mat_b(i) = -1 / (D%z(i+1)-D%z(i)) / (D%z_u(i) - D%z_u(i-1))
        D%mat_b(i) = D%mat_b(i) - 1 / (D%z(i)-D%z(i-1)) / (D%z_u(i) - D%z_u(i-1))   
      end do
      if (D%gnz>1) then
        i = D%gnz
        D%mat_b(D%gnz) = -1 / (D%z(i)-D%z(i-1)) / (D%z_u(i) - D%z_u(i-1))
      end if

      !Dirichlet: the flux to the zero value at the face z_u(0) or z_u(n)
      n = D%gnz
      if (D%BCs(5)==PoisFFT_DirichletStag) &
        D%mat_b(1) = D%mat_b(1) - 1 / (D%z(1)-D%z_u(0)) / (D%z_u(1) - D%z_u(0))
      if (D%BCs(6)==PoisFFT_DirichletStag) &
        D%mat_b(n) = D%mat_b(n) - 1 / (D%z_u(n)-D%z(n)) / (D%z_u(n) - D%z_u(n-1))
      
      do i = 1, D%gnz-1
        D%mat_c(i) = 1 / (D%z(i+1)-D%z(i)) / (D%z_u(i) - D%z_u(i-1))
      end do

    contains

      !returns the error in ierr and errmsg if ierr is present, stops otherwise
      subroutine report_error(code, detail)
        integer, intent(in) :: code
        integer, intent(in), optional :: detail
        character(:), allocatable :: msg

        select case (code)
        case (PoisFFT_ErrorBCs)
          msg = "Unknown combination of boundary conditions."
        case default
          if (detail==1) then
            msg = "The arrays z(1:) and z_u(0:) must cover the whole grid in z."
          else
            msg = "The coordinates z and z_u must be increasing."
          end if
        end select

        if (present(ierr)) then
          ierr = code
          if (present(errmsg)) errmsg = msg
        else
          !literal messages, a variable stop code requires Fortran 2018
          select case (code)
          case (PoisFFT_ErrorBCs)
            stop "Unknown combination of boundary conditions."
          case default
            if (detail==1) then
              stop "The arrays z(1:) and z_u(0:) must cover the whole grid in z."
            else
              stop "The coordinates z and z_u must be increasing."
            end if
          end select
        end if
      end subroutine

    end subroutine PoisFFT_Solver3D_nonuniform_z_Init
    
    


    subroutine PoisFFT_Solver3D_nonuniform_z__Execute(D,Phi,RHS)
      type(PoisFFT_Solver3D_nonuniform_z), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)

      integer   :: ngPhi(3), ngRHS(3)

      ngPhi = (ubound(Phi)-[D%nx,D%ny,D%nz])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny,D%nz])/2

      if (.not.all(D%BCs(1:4)==PoisFFT_Periodic)) then

        call PoisFFT_Solver3D_nonuniform_z_RealXY(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else

        call PoisFFT_Solver3D_nonuniform_z_PeriodicXY(D,&
                Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                    ngPhi(2)+1:ngPhi(2)+D%ny,&
                    ngPhi(3)+1:ngPhi(3)+D%nz),&
                RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                    ngRHS(2)+1:ngRHS(2)+D%ny,&
                    ngRHS(3)+1:ngRHS(3)+D%nz))

      endif

    end subroutine PoisFFT_Solver3D_nonuniform_z__Execute
    
    
    
    
    subroutine PoisFFT_Solver3D_nonuniform_z__Finalize(D)
      type(PoisFFT_Solver3D_nonuniform_z), intent(inout) :: D
      integer :: i

      call Finalize(D%PoisFFT_Solver3D)

      !not all of them are allocated when the constructor returned an error
      if (allocated(D%z)) deallocate(D%z)
      if (allocated(D%z_u)) deallocate(D%z_u)

      if (allocated(D%mat_a)) deallocate(D%mat_a)
      if (allocated(D%mat_b)) deallocate(D%mat_b)
      if (allocated(D%mat_c)) deallocate(D%mat_c)

    endsubroutine PoisFFT_Solver3D_nonuniform_z__Finalize
    
    



    
    
    
    
    
    
    
    
    
#include "poisfft-solvers-inc.f90"

#include "poisfft-solvers-nonuniform_z-inc.f90"






    function eigenvalues(f, BCs, L, n, gn, off) result(res)
      procedure(eig_fn_spectral) :: f
      integer, intent(in) :: BCs(2)
      real(RP), intent(in) :: L
      integer, intent(in) :: n, gn, off
      real(RP) :: res(n)
      integer :: i
      real(RP) :: dkx, dkx_h
      
      dkx = pi * grid_dx(gn, L, BCs) / L
      dkx_h = dkx / 2
      
      if (all(BCs==PoisFFT_Periodic)) then
        !the wave number of the DFT mode m = i-1+off is min(m, gn-m)
        do i = 1, n
          if (i-1+off<=gn/2) then
            res(i) = f((i-1+off)*dkx)
          else
            res(i) = f((gn-(i-1+off))*dkx)
          end if
        end do
      else if (all(BCs==PoisFFT_Dirichlet)) then
        forall(i=1:n) res(i) = f((i+off)*dkx_h)
      else if (all(BCs==PoisFFT_DirichletStag)) then
        forall(i=1:n) res(i) = f((i+off)*dkx_h)
      else if (all(BCs==PoisFFT_Neumann)) then
        forall(i=1:n) res(i) = f((i-1+off)*dkx_h)
      else if (all(BCs==PoisFFT_NeumannStag)) then
        forall(i=1:n) res(i) = f((i-1+off)*dkx_h)
      else if (all(BCs==[PoisFFT_NeumannStag,PoisFFT_DirichletStag])) then
        forall(i=1:n) res(i) = f((i+off-0.5_RP)*dkx_h)
      else if (all(BCs==[PoisFFT_DirichletStag,PoisFFT_NeumannStag])) then
        forall(i=1:n) res(i) = f((i+off-0.5_RP)*dkx_h)
      else if (all(BCs==[PoisFFT_Neumann,PoisFFT_Dirichlet])) then
        forall(i=1:n) res(i) = f((i+off-0.5_RP)*dkx_h)
      else if (all(BCs==[PoisFFT_Dirichlet,PoisFFT_Neumann])) then
        forall(i=1:n) res(i) = f((i+off-0.5_RP)*dkx_h)
      end if

      res = res / (grid_dx(gn, L, BCs))**2
      
    end function


    pure real(RP) function eig_fn_spectral(x) result(f)
      real(RP), intent(in) :: x
      f = -(2 * x)**2
    end function
  
    pure real(RP) function eig_fn_FD2(x) result(f)
      real(RP), intent(in) :: x
      f = -(2 * sin(x))**2
    end function
  
    pure real(RP) function eig_fn_FD4(x) result(f)
      real(RP), intent(in) :: x
      f = -((sin(3*x) - 27*sin(x)) / 12)**2
    end function
    
      
    function grid_dx(gn, L, BCs) result(res)
      real(RP) :: res
      integer, intent(in) :: gn
      real(RP), intent(in) :: L
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Periodic)) then
        res = L / gn
      else if (all(BCs==PoisFFT_Dirichlet)) then
        res = L / (gn+1)
      else if (all(BCs==PoisFFT_DirichletStag .or. &
                   BCs==PoisFFT_NeumannStag)) then
        res = L / gn
      else if (all(BCs==PoisFFT_Neumann)) then
        res = L / (gn-1)
      else if (all(BCs==PoisFFT_Dirichlet .or. &
                   BCs==PoisFFT_Neumann)) then
        res = L / gn
      end if
    end function

    function norm_factor(gn,BCs) result(res)
      real(RP) :: res
      integer, intent(in) :: gn
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Periodic)) then
        res = gn
      else if (all(BCs==PoisFFT_Dirichlet)) then
        res = 2 * (gn+1)
      else if (all(BCs==PoisFFT_Neumann)) then
        res = 2 * (gn-1)
      else if (all(BCs==PoisFFT_DirichletStag .or. &
                   BCs==PoisFFT_NeumannStag)) then
        res = 2 * gn
      else if (all(BCs==PoisFFT_Dirichlet .or. &
                   BCs==PoisFFT_Neumann)) then
        res = 2 * gn
      end if
    end function
    
    function real_transform_type_forward(BCs) result(res)
      integer :: res
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Periodic)) then
          !the halfcomplex order is compatible with the order of eigenvalues()
          res = FFT_RealToHalfComplex
      else if (all(BCs==PoisFFT_Dirichlet)) then
          res = FFT_RealOdd00
      else if (all(BCs==PoisFFT_DirichletStag)) then
          res = FFT_RealOdd10
      else if (all(BCs==PoisFFT_Neumann)) then
          res = FFT_RealEven00
      else if (all(BCs==PoisFFT_NeumannStag)) then
          res = FFT_RealEven10
      else if (all(BCs==[PoisFFT_NeumannStag,PoisFFT_DirichletStag])) then
          res = FFT_RealEven11
      else if (all(BCs==[PoisFFT_DirichletStag,PoisFFT_NeumannStag])) then
          res = FFT_RealOdd11
      else if (all(BCs==[PoisFFT_Neumann,PoisFFT_Dirichlet])) then
          res = FFT_RealEven01
      else if (all(BCs==[PoisFFT_Dirichlet,PoisFFT_Neumann])) then
          res = FFT_RealOdd01
      else
          res = FFT_Unsupported
      end if
    end function
    
    function real_transform_type_backward(BCs) result(res)
      integer :: res
      integer, intent(in) :: BCs(2)
      
      if (all(BCs==PoisFFT_Periodic)) then
          res = FFT_HalfComplexToReal
      else if (all(BCs==PoisFFT_Dirichlet)) then
          res = FFT_RealOdd00
      else if (all(BCs==PoisFFT_DirichletStag)) then
          res = FFT_RealOdd01
      else if (all(BCs==PoisFFT_Neumann)) then
          res = FFT_RealEven00
      else if (all(BCs==PoisFFT_NeumannStag)) then
          res = FFT_RealEven01
      else if (all(BCs==[PoisFFT_NeumannStag,PoisFFT_DirichletStag])) then
          res = FFT_RealEven11
      else if (all(BCs==[PoisFFT_DirichletStag,PoisFFT_NeumannStag])) then
          res = FFT_RealOdd11
      else if (all(BCs==[PoisFFT_Neumann,PoisFFT_Dirichlet])) then
          res = FFT_RealEven10
      else if (all(BCs==[PoisFFT_Dirichlet,PoisFFT_Neumann])) then
          res = FFT_RealOdd10
      else
          res = FFT_Unsupported
      end if
    end function


    !Stops for a pair of boundary conditions in any direction that has no
    !transform (e.g. periodic with non-periodic, or a regular with a staggered
    !Dirichlet condition). The supported combinations of the directions are
    !checked by the Init procedures themselves.
    subroutine check_BC_pairs(BCs)
      integer, intent(in) :: BCs(:)

      if (.not.BC_pairs_supported(BCs)) &
        stop "Unknown combination of boundary conditions."
    end subroutine

    logical function BC_pairs_supported(BCs) result(res)
      integer, intent(in) :: BCs(:)
      integer :: i

      res = .true.
      do i = 1, size(BCs), 2
        if (real_transform_type_forward(BCs(i:i+1))==FFT_Unsupported) res = .false.
      end do
    end function

    

    
#ifdef MPI     
    subroutine Init_MPI_Buffers(D, dir, complex)
      
      class(PoisFFT_Solver3D), intent(inout), target :: D
      integer, intent(in) :: dir
      logical, intent(in), optional :: complex
      type(mpi_vars_1d), pointer :: mpi
      integer :: i, ie
      logical :: complex_buffers

      complex_buffers = .false.
      if (present(complex)) complex_buffers = complex

      mpi => D%Solvers1D(dir)%mpi
    
      allocate(mpi%snxs(mpi%np))
      allocate(mpi%snzs(mpi%np))

      allocate(mpi%rnxs(mpi%np))
      allocate(mpi%rnzs(mpi%np))

      allocate(mpi%sumrnzs(mpi%np))

      allocate(mpi%sdispls(mpi%np))
      allocate(mpi%scounts(mpi%np))

      allocate(mpi%rdispls(mpi%np))
      allocate(mpi%rcounts(mpi%np))
        
      if (dir==2) then
        mpi%snxs(1:mpi%np-1) = (D%nx / mpi%np)
        mpi%snxs(mpi%np) = D%nx - sum(mpi%snxs(1:mpi%np-1))
        mpi%snzs = D%ny
        mpi%scounts = mpi%snxs*D%ny*D%nz
        mpi%sdispls(1) = 0
        do i = 2, mpi%np
          mpi%sdispls(i) = mpi%sdispls(i-1) + mpi%scounts(i-1)
        end do

        call MPI_AllToAll(mpi%snxs, 1, MPI_INTEGER, &
                          mpi%rnxs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        if (.not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))) &
          stop "PoisFFT internal error: .not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))"
        call MPI_AllToAll(mpi%snzs, 1, MPI_INTEGER, &
                          mpi%rnzs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        call MPI_AllToAll(mpi%scounts, 1, MPI_INTEGER, &
                          mpi%rcounts, 1, MPI_INTEGER, &
                          mpi%comm, ie)

        mpi%rdispls(1) = 0
        do i = 2, mpi%np
          mpi%rdispls(i) = mpi%rdispls(i-1) + mpi%rcounts(i-1)
        end do

        do i = 1, mpi%np
          mpi%sumrnzs(i) = sum(mpi%rnzs(1:i-1))
        end do

        if (complex_buffers) then
          allocate(mpi%ctmp1(1:D%ny,1:D%nz,1:D%nx))
          allocate(mpi%ctmp2(0:sum(mpi%rcounts)-1))
          allocate(mpi%cwork(sum(mpi%rnzs), D%nz, mpi%rnxs(1)))
        else
          allocate(mpi%tmp1(1:D%ny,1:D%nz,1:D%nx))
          allocate(mpi%tmp2(0:sum(mpi%rcounts)-1))
          allocate(mpi%rwork(sum(mpi%rnzs), D%nz, mpi%rnxs(1)))
        end if
        
      else if (dir==3) then
      
        mpi%snxs(1:mpi%np-1) = (D%nx / mpi%np)
        mpi%snxs(mpi%np) = D%nx - sum(mpi%snxs(1:mpi%np-1))
        mpi%snzs = D%nz
        mpi%scounts = mpi%snxs*D%ny*D%nz
        mpi%sdispls(1) = 0
        do i = 2, mpi%np
          mpi%sdispls(i) = mpi%sdispls(i-1) + mpi%scounts(i-1)
        end do

        call MPI_AllToAll(mpi%snxs, 1, MPI_INTEGER, &
                          mpi%rnxs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        if (.not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))) &
          stop "PoisFFT internal error: .not.all(mpi%rnxs(2:mpi%np)==mpi%rnxs(1))"
        call MPI_AllToAll(mpi%snzs, 1, MPI_INTEGER, &
                          mpi%rnzs, 1, MPI_INTEGER, &
                          mpi%comm, ie)
        call MPI_AllToAll(mpi%scounts, 1, MPI_INTEGER, &
                          mpi%rcounts, 1, MPI_INTEGER, &
                          mpi%comm, ie)

        mpi%rdispls(1) = 0
        do i = 2, mpi%np
          mpi%rdispls(i) = mpi%rdispls(i-1) + mpi%rcounts(i-1)
        end do

        do i = 1, mpi%np
          mpi%sumrnzs(i) = sum(mpi%rnzs(1:i-1))
        end do

        if (complex_buffers) then
          allocate(mpi%ctmp1(1:D%nz,1:D%ny,1:D%nx))
          allocate(mpi%ctmp2(0:sum(mpi%rcounts)-1))
          allocate(mpi%cwork(sum(mpi%rnzs), D%ny, mpi%rnxs(1)))
        else
          allocate(mpi%tmp1(1:D%nz,1:D%ny,1:D%nx))
          allocate(mpi%tmp2(0:sum(mpi%rcounts)-1))
          allocate(mpi%rwork(sum(mpi%rnzs), D%ny, mpi%rnxs(1)))
        end if
      else
        stop "Not implemented."
      end if
    end subroutine Init_MPI_Buffers


    !Initialization of PoisFFT_Solver3D_Separable_MPI, used for the
    !combinations of boundary conditions without a specialized MPI solver.
    !Every direction uses a real 1D transform (halfcomplex for periodic),
    !Solvers1D(dir + 3*(thread-1)) as in the MPI version of PNsNs.
    subroutine PoisFFT_Solver3D_Separable_MPI_Init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: dir, i

      allocate(D%Solvers1D(3*D%nthreads))

      do dir = 1, 3
        D%Solvers1D(dir) = PoisFFT_Solver1D_From3D(D,dir)

        call allocate_fftw_real(D%Solvers1D(dir))

        D%Solvers1D(dir)%forward = PoisFFT_Plan1D(D%Solvers1D(dir), &
                                     [real_transform_type_forward(D%BCs(2*dir-1:2*dir))])
        D%Solvers1D(dir)%backward = PoisFFT_Plan1D(D%Solvers1D(dir), &
                                     [real_transform_type_backward(D%BCs(2*dir-1:2*dir))])

        do i = dir + 3, dir + 3*(D%nthreads-1), 3
          D%Solvers1D(i) = D%Solvers1D(dir)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call allocate_fftw_real(D%Solvers1D(i))
        end do
      end do

      call Init_MPI_Buffers(D, 2)
      call Init_MPI_Buffers(D, 3)
    end subroutine PoisFFT_Solver3D_Separable_MPI_Init
#endif

#undef RP
#undef CP
#undef _RP
#undef _CP
#undef MPI_RP
#undef MPI_CP
