#if (PREC==2)

#define RP DRP
#define CP DCP
#define MPI_RP MPI_DOUBLE_PRECISION

#define _RP _DRP

#else

#define RP SRP
#define CP SCP
#define MPI_RP MPI_REAL

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
            Finalize, Execute

  real(RP), parameter, private :: pi = 3.141592653589793238462_RP


  interface PoisFFT_Solver1D
    module procedure PoisFFT_Solver1D__New
  end interface

  interface PoisFFT_Solver2D
    module procedure PoisFFT_Solver2D__New
  end interface

  interface PoisFFT_Solver3D
    module procedure PoisFFT_Solver3D__New
  end interface

  interface Finalize
    module procedure PoisFFT_Solver1D__Finalize
    module procedure PoisFFT_Solver2D__Finalize
    module procedure PoisFFT_Solver3D__Finalize
  end interface Finalize


  interface Init
    module procedure PoisFFT_Solver1D_Init
    module procedure PoisFFT_Solver2D_Init
    module procedure PoisFFT_Solver3D_Init
  end interface Init

  interface Execute
    module procedure PoisFFT_Solver1D__Execute
    module procedure PoisFFT_Solver2D__Execute
    module procedure PoisFFT_Solver3D__Execute
  end interface Execute


  contains


    function PoisFFT_Solver3D__New(nxyz,Lxyz,BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver3D) :: D

      integer, intent(in)   :: nxyz(3)
      real(RP), intent(in)  :: Lxyz(3)
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(3)
      integer, intent(in), optional :: offs(3)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads

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
      
      if (D%BCs(1)/=D%BCs(2) .or. &
          D%BCs(3)/=D%BCs(4) .or. &
          D%BCs(5)/=D%BCs(6)) then
        stop "Both boundary consitions in one direction must be identical."
      end if
      
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

      !create fftw plans and allocate working arrays
      call Init(D)
    end function PoisFFT_Solver3D__New
    
    


    subroutine PoisFFT_Solver3D_Init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: real_forw, real_back, i
      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel

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
        
        real_forw = real_transform_type_forward(D%BCs(1))
        real_back = real_transform_type_backward(D%BCs(1))

        D%forward = PoisFFT_Plan3D(D, [(real_forw, i=1,3)])
        D%backward = PoisFFT_Plan3D(D, [(real_back, i=1,3)])

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. &
               (all(D%BCs(5:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(5:6)==PoisFFT_NeumannStag) )) then

#ifdef MPI
        allocate(D%Solvers1D(3:2+D%nthreads))

        D%Solvers1D(3) = PoisFFT_Solver1D_From3D(D,3)

        call allocate_fftw_real(D%Solvers1D(3))

        D%Solvers1D(3)%forward = PoisFFT_Plan1D(D%Solvers1D(3), [FFT_RealEven10])
        D%Solvers1D(3)%backward = PoisFFT_Plan1D(D%Solvers1D(3), [FFT_RealEven01])

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

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_RealEven10])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_RealEven01])

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

      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               all(D%BCs(5:6)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:4)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:4)==PoisFFT_NeumannStag) )) then

        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,2)

        call allocate_fftw_real(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_RealEven10])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D(D%Solvers1D(1), [FFT_RealEven01])

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
        
      else if (all(D%BCs(1:2)==PoisFFT_Periodic) .and. &
               (all(D%BCs(3:6)==PoisFFT_Neumann) .or. &
                all(D%BCs(3:6)==PoisFFT_NeumannStag) )) then

        real_forw = real_transform_type_forward(D%BCs(3))
        real_back = real_transform_type_backward(D%BCs(3))
        
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
      else
        stop "Unknown combination of boundary conditions."
      endif


      if (D%approximation==2) then
        D%denomx = eigenvalues_FD2(D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues_FD2(D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues_FD2(D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      else
        D%denomx = eigenvalues_spectral(D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues_spectral(D%BCs(3:4), D%Ly, D%ny, D%gny, D%offy)
        D%denomz = eigenvalues_spectral(D%BCs(5:6), D%Lz, D%nz, D%gnz, D%offz)
      end if
      
    contains
    
      function real_transform_type_forward(bc) result(res)
        integer :: res
        integer, intent(in) :: bc
        
        select case (bc)
          case (PoisFFT_Dirichlet)
            res = FFT_RealOdd00
          case (PoisFFT_DirichletStag)
            res = FFT_RealOdd10
          case (PoisFFT_Neumann)
            res = FFT_RealEven00
          case (PoisFFT_NeumannStag)
            res = FFT_RealEven10
          case default
            res = -1
        end select
      end function
      
      function real_transform_type_backward(bc) result(res)
        integer :: res
        integer, intent(in) :: bc
        
        select case (bc)
          case (PoisFFT_Dirichlet)
            res = FFT_RealOdd00
          case (PoisFFT_DirichletStag)
            res = FFT_RealOdd01
          case (PoisFFT_Neumann)
            res = FFT_RealEven00
          case (PoisFFT_NeumannStag)
            res = FFT_RealEven01
          case default
            res = -1
        end select
      end function
      
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
      type(PoisFFT_Solver3D), intent(in) :: D3D
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
    end function PoisFFT_Solver1D_From3D



    function PoisFFT_Solver2D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver2D)             :: D
      type(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction
#ifdef MPI
      integer :: ie, dims
      
      call MPI_Cartdim_get(D3D%mpi%comm, dims, ie)
      if (ie/=0) stop "Error executing MPI_Cartdim_get."
      
      if (dims==2) then
        !We see the dimensions reversed in Fortran!
        if (direction==1) then
          D%mpi%comm = D3D%mpi%comm
          D%mpi%comm_dim = 2
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
    end function PoisFFT_Solver2D_From3D




    function PoisFFT_Solver2D__New(nxyz,Lxyz,BCs,approximation, &
                                  gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver2D) :: D

      integer, intent(in)   :: nxyz(2)
      real(RP), intent(in)  :: Lxyz(2)
      integer, intent(in)   :: bcs(4)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(2)
      integer, intent(in), optional :: offs(2)
      integer, intent(in), optional  :: mpi_comm
      integer, intent(in), optional :: nthreads

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

      !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver2D__New




    subroutine Poisfft_Solver2D_Init(D)
      type(PoisFFT_Solver2D), intent(inout) :: D
      integer :: i
        
      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2)) * &
                      norm_factor(D%gny,D%BCs(3:4))
      
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call allocate_fftw_complex(D)

        D%forward = PoisFFT_Plan2D(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan2D(D, [FFT_Complex, FFTW_BACKWARD])

      else if (all(D%BCs==PoisFFT_Dirichlet)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan2D(D, [(FFT_RealOdd00, i=1,2)])
        D%backward = PoisFFT_Plan2D(D, [(FFT_RealOdd00, i=1,2)])

       else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan2D(D, [(FFT_RealOdd10, i=1,2)])
        D%backward = PoisFFT_Plan2D(D, [(FFT_RealOdd01, i=1,2)])

      else if (all(D%BCs==PoisFFT_Neumann)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan2D(D, [(FFT_RealEven00, i=1,2)])

        D%backward = PoisFFT_Plan2D(D, [(FFT_RealEven00, i=1,2)])

      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan2D(D, [(FFT_RealEven10, i=1,2)])
        D%backward = PoisFFT_Plan2D(D, [(FFT_RealEven01, i=1,2)])

      endif
      
      if (D%approximation==2) then
        D%denomx = eigenvalues_FD2(D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues_FD2(D%BCs(2:4), D%Ly, D%ny, D%gny, D%offy)
      else
        D%denomx = eigenvalues_spectral(D%BCs(1:2), D%Lx, D%nx, D%gnx, D%offx)
        D%denomy = eigenvalues_spectral(D%BCs(2:4), D%Ly, D%ny, D%gny, D%offy)
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

      else if (all(D%BCs==PoisFFT_Dirichlet) .or. &
               all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver2D_FullDirichlet(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny))

      else if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver2D_FullNeumann(D,&
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
                                  gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver1D) :: D

      integer, intent(in)   :: nxyz(1)
      real(RP), intent(in)  :: Lxyz(1)
      integer, intent(in)   :: bcs(2)
      integer, intent(in), optional :: approximation
      integer, intent(in), optional :: gnxyz(1)
      integer, intent(in), optional :: offs(1)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads
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

       !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver1D__New




    subroutine Poisfft_Solver1D_Init(D)
      type(PoisFFT_Solver1D), intent(inout) :: D
      
      D%norm_factor = norm_factor(D%gnx,D%BCs(1:2))
      
      allocate(D%denomx(D%gnx))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call allocate_fftw_complex(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan1D(D, [FFT_Complex, FFTW_BACKWARD])

      else if (all(D%BCs==PoisFFT_Dirichlet)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealOdd00])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealOdd00])
        
      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealOdd10])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealOdd01])
        
      else if (all(D%BCs==PoisFFT_Neumann)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealEven00])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealEven00])
        
      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call allocate_fftw_real(D)

        D%forward = PoisFFT_Plan1D(D, [FFT_RealEven10])
        D%backward = PoisFFT_Plan1D(D, [FFT_RealEven01])
        
      endif
      
      if (D%approximation==2) then
        D%denomx = eigenvalues_FD2(D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
      else
        D%denomx = eigenvalues_spectral(D%BCs(1:2), D%Lx, D%gnx, D%gnx, D%offx)
      end if

    end subroutine PoisFFT_Solver1D_Init
    
    
    
    
    subroutine PoisFFT_Solver1D__Execute(D,Phi,RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)

      integer   :: ngPhi(1), ngRHS(1)


      ngPhi = (ubound(Phi)-D%nx)/2
      ngRHS = (ubound(RHS)-D%nx)/2


      if (all(D%BCs==PoisFFT_Periodic)) then

        call PoisFFT_Solver1D_FullPeriodic(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))

      else if (all(D%BCs==PoisFFT_Dirichlet) .or. &
               all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver1D_FullDirichlet(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))

      else if (all(D%BCs==PoisFFT_Neumann) .or. &
               all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver1D_FullNeumann(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))
                 
      endif

    end subroutine PoisFFT_Solver1D__Execute




    subroutine PoisFFT_Solver1D__Finalize(D)
      type(PoisFFT_Solver1D), intent(inout) :: D

      call Finalize(D%forward)
      call Finalize(D%backward)

      if (associated(D%rwork)) call deallocate_fftw(D%rwork)
      if (associated(D%cwork)) call deallocate_fftw(D%cwork)

    endsubroutine PoisFFT_Solver1D__Finalize





#include "poisfft-solvers-inc.f90"





    function eigenvalues_spectral(BCs, L, n, gn, off) result(res)
      integer, intent(in) :: BCs(2)
      real(RP), intent(in) :: L
      integer, intent(in) :: n, gn, off
      real(RP) :: res(n)
      integer :: i
      real(RP) :: dkx, dkx_h
      
      dkx_h = pi/(gn)
      dkx = 2 * dkx_h
      
      if (all(BCs==PoisFFT_Periodic)) then
        do i = 1, n
          if (i+off<gn/2) then
            res(i) = -4*pi**2*(i-1+off)**2
          else
            res(i) = -4*pi**2*(gn-i+1+off)**2
          end if
        end do
      else if (all(BCs==PoisFFT_Dirichlet)) then
        forall(i=1:n) res(i) = -pi**2*(i+off)**2
      else if (all(BCs==PoisFFT_DirichletStag)) then
        forall(i=1:n) res(i) = -pi**2*(i+off)**2
      else if (all(BCs==PoisFFT_Neumann)) then
        forall(i=1:n) res(i) = -pi**2*(i-1+off)**2
      else if (all(BCs==PoisFFT_NeumannStag)) then
        forall(i=1:n) res(i) = -pi**2*(i-1+off)**2
      endif
      
      res = res / L**2
    end function

    
    

    function eigenvalues_FD2(BCs, L, n, gn, off) result(res)
      integer, intent(in) :: BCs(2)
      real(RP), intent(in) :: L
      integer, intent(in) :: n, gn, off
      real(RP) :: res(n)
      integer :: i
      real(RP) :: dkx, dkx_h
      
      dkx_h = pi/(gn)
      dkx = 2 * dkx_h
      
      if (all(BCs==PoisFFT_Periodic)) then
        forall(i=1:n) res(i) = 2 * (cos((i-1+off)*dkx)-1.0_RP)
      else if (all(BCs==PoisFFT_Dirichlet)) then
        forall(i=1:n) res(i) = 2 * (cos((i+off)*pi/(gn+1))-1.0_RP)
      else if (all(BCs==PoisFFT_DirichletStag)) then
        forall(i=1:n) res(i) = 2 * (cos((i+off)*dkx_h)-1.0_RP)
      else if (all(BCs==PoisFFT_Neumann)) then
        forall(i=1:n) res(i) = 2 * (cos((i-1+off)*pi/(gn-1))-1.0_RP)
      else if (all(BCs==PoisFFT_NeumannStag)) then
        forall(i=1:n) res(i) = 2 * (cos((i-1+off)*dkx_h)-1.0_RP)
      end if

      res = res / (grid_dx(gn, L, BCs))**2
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
      else if (all(BCs==PoisFFT_DirichletStag)) then
        res = L / gn
      else if (all(BCs==PoisFFT_Neumann)) then
        res = L / (gn-1)
      else if (all(BCs==PoisFFT_NeumannStag)) then
        res = L / gn
      else
        res = 0
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
      else if (all(BCs==PoisFFT_DirichletStag)) then
        res = 2 * gn
      else if (all(BCs==PoisFFT_Neumann)) then
        res = 2 * (gn-1)
      else if (all(BCs==PoisFFT_NeumannStag)) then
        res = 2 *  gn
      else
        res = 0
      end if
    end function
    
#ifdef MPI     
    subroutine Init_MPI_Buffers(D, dir)
      interface
        subroutine MPI_ALLTOALL(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,&
                                RECVTYPE, COMM, IERROR)
            INTEGER   SENDBUF(*), RECVBUF(*)
            INTEGER   SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
            INTEGER   COMM, IERROR
        end subroutine
      end interface
      
      type(PoisFFT_Solver3D), intent(inout), target :: D
      integer, intent(in) :: dir
      type(mpi_vars_1d), pointer :: mpi
      integer :: i, ie
      
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

        allocate(mpi%tmp1(1:D%ny,1:D%nz,1:D%nx))
        allocate(mpi%tmp2(0:sum(mpi%rcounts)-1))
        allocate(mpi%rwork(sum(mpi%rnzs), D%nz, mpi%rnxs(1)))
        
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

        allocate(mpi%tmp1(1:D%nz,1:D%ny,1:D%nx))
        allocate(mpi%tmp2(0:sum(mpi%rcounts)-1))
        allocate(mpi%rwork(sum(mpi%rnzs), D%ny, mpi%rnxs(1)))
      else
        stop "Not implemented."
      end if
    end subroutine Init_MPI_Buffers
#endif

#undef RP
#undef CP
#undef _RP
#undef _CP
#undef MPI_RP
