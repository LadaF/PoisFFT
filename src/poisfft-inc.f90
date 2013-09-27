  use iso_c_binding
  !$ use omp_lib
  use PoisFFT_Precisions
  use PoisFFT_Parameters

  implicit none

  private
  public :: PoisFFT_Solver1D, PoisFFT_Solver2D, PoisFFT_Solver3D, &
            DeallocateData, Execute, New

  real(RP), parameter, private :: pi = 3.141592653589793238462_RP


  interface New
    module procedure PoisFFT_Solver1D_New
    module procedure PoisFFT_Solver2D_New
    module procedure PoisFFT_Solver3D_New
  end interface New

  interface DeallocateData
    module procedure PoisFFT_Solver1D_DeallocateData
    module procedure PoisFFT_Solver2D_DeallocateData
    module procedure PoisFFT_Solver3D_DeallocateData
  end interface DeallocateData


  interface Init
    module procedure PoisFFT_Solver1D_Init
    module procedure PoisFFT_Solver2D_Init
    module procedure PoisFFT_Solver3D_Init
  end interface Init

  interface Execute
    module procedure PoisFFT_Solver1D_Execute
    module procedure PoisFFT_Solver2D_Execute
    module procedure PoisFFT_Solver3D_Execute
  end interface Execute


  contains

    subroutine PoisFFT_Solver3D_Execute(D,Phi,RHS)
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

      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver3D_FullDirichletStag(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver3D_FullNeumannStag(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

      else if (all(D%BCs(1:4)==PoisFFT_Periodic).and.all(D%BCs(5:6)==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver3D_PPNs(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

       endif

    end subroutine PoisFFT_Solver3D_Execute



    subroutine PoisFFT_Solver3D_New(D,nx,ny,nz,dx,dy,dz,BCs,gnxyz,offs,mpi_comm,nthreads)
#ifdef MPI
  use mpi
#endif
      type(PoisFFT_Solver3D), intent(out) :: D

      integer, intent(in)   :: nx, ny, nz
      real(RP), intent(in)  :: dx, dy, dz
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: nthreads
      integer(c_size_t), intent(in), optional :: gnxyz(3)
      integer(c_size_t), intent(in), optional :: offs(3)
      integer(c_int32_t), intent(in), optional  :: mpi_comm

      D%dx = dx
      D%dy = dy
      D%dz = dz

      D%nx = nx
      D%ny = ny
      D%nz = nz
      
      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
        D%gny = gnxyz(2)
        D%gnz = gnxyz(3)
      else
        D%gnx = nx
        D%gny = ny
        D%gnz = nz
      end if

      if (present(offs)) then
        D%offx = offs(1)
        D%offy = offs(2)
        D%offz = offs(3)
      end if

      D%cnt = D%nx * D%ny * D%nz

      D%BCs = BCs
      
#ifdef MPI      
      if (present(mpi_comm)) then
        D%mpi_comm = mpi_comm
      else
        stop "No PFFT comm present in PoisFFT_Solver3D_New."
      end if
#endif

      if (present(nthreads)) then
        D%nthreads = nthreads
      else
        D%nthreads = 1
      endif

      !create fftw plans and allocate working arrays
      call Init(D)
    end subroutine PoisFFT_Solver3D_New


    function PoisFFT_Solver1D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver1D)             :: D
      type(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction

      if (direction==1) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%BCs = D3D%BCs(1:2)
      else if (direction==2) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%BCs = D3D%BCs(3:4)
      else
        D%dx = D3D%dz
        D%nx = D3D%nz
        D%gnx = D3D%gnz
        D%BCs = D3D%BCs(5:6)
      endif

      D%cnt = D%nx
    end function PoisFFT_Solver1D_From3D


    function PoisFFT_Solver2D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver2D)             :: D
      type(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction

      if (direction==1) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%BCs = D3D%BCs(3:6)
      else if (direction==2) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%BCs = D3D%BCs([1,2,5,6])
      else
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%dy = D3D%dy
        D%ny = D3D%ny
        D%gny = D3D%gny
        D%BCs = D3D%BCs(1:4)
      endif

      D%cnt = D%nx * D%ny
    end function PoisFFT_Solver2D_From3D


    subroutine PoisFFT_Solver1D_DeallocateData(D)
      type(PoisFFT_Solver1D), intent(inout) :: D

      call Destroy(D%forward)
      call Destroy(D%backward)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

    endsubroutine PoisFFT_Solver1D_DeallocateData



    subroutine PoisFFT_Solver2D_DeallocateData(D)
      type(PoisFFT_Solver2D), intent(inout) :: D

      call Destroy(D%forward)
      call Destroy(D%backward)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

    endsubroutine PoisFFT_Solver2D_DeallocateData



    subroutine PoisFFT_Solver3D_DeallocateData(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: i

      call Destroy(D%forward)
      call Destroy(D%backward)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

      if (allocated(D%Solvers1D)) then
        do i=1,size(D%Solvers1D)
          call DeallocateData(D%Solvers1D(i))
        end do
      endif

      if (allocated(D%Solvers2D)) then
        do i=1,size(D%Solvers2D)
          call DeallocateData(D%Solvers2D(i))
        end do
      endif

      if (allocated(D%Solvers2D)) then
        do i=1,size(D%Solvers2D)
          call DeallocateData(D%Solvers2D(i))
        end do
      endif
    endsubroutine PoisFFT_Solver3D_DeallocateData



    subroutine PoisFFT_Solver3D_FullPeriodic(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      real(RP) :: dx2, dy2, dz2
      real(RP) :: dkx, dky, dkz
      real(RP) :: denomx(D%nx), denomy(D%ny), denomz(D%nz)
      integer i,j,k

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2
      dz2 = 1._RP/D%dz**2

      dkx = 2._RP*pi/(D%gnx)
      dky = 2._RP*pi/(D%gny)
      dkz = 2._RP*pi/(D%gnz)

      !$omp parallel private(i,j,k)

      !$omp workshare
      D%cwork = cmplx(RHS,0._RP,CP)
      !$omp end workshare

      !$omp end parallel
#ifdef MPI
      call PoisFFT_Plan_Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%cwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp sections
      !$omp section
      forall(i=1:D%nx)  denomx(i) = (cos((i-1+D%offx)*dkx)-1.0_RP) * dx2 * 2
      !$omp section
      forall(j=1:D%ny)  denomy(j) = (cos((j-1+D%offy)*dky)-1.0_RP) * dy2 * 2
      !$omp section
      forall(k=1:D%nz)  denomz(k) = (cos((k-1+D%offz)*dkz)-1.0_RP) * dz2 * 2
      !$omp end sections

      if (D%offx==0.and.D%offy==0.and.D%offz==0) then
#define xwork cwork
#include "loop_nest_3d.f90"
#undef xwork
        !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
        ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
      
        !$omp single
        D%cwork(1,1,1) = 0
        !$omp end single
      else
        !$omp do
        do k=1,D%nz
          do j=1,D%ny
            do i=1,D%nx
              D%cwork(i,j,k) = D%cwork(i,j,k) / (denomx(i) + denomy(j) + denomz(k))
            end do
          end do
        end do
        !$omp end do
      end if

      !$omp end parallel
#ifdef MPI
      call PoisFFT_Plan_Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%cwork)
#endif
      !$omp parallel

      !$omp workshare
       Phi = real(D%cwork,RP) / (D%gnx*D%gny*D%gnz)
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullPeriodic




    subroutine PoisFFT_Solver3D_FullDirichletStag(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      real(RP) :: dx2, dy2, dz2
      real(RP) :: dkx_h, dky_h, dkz_h !half of ussual dkx
      real(RP) :: denomx(D%nx), denomy(D%ny), denomz(D%nz)
      integer i,j,k

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2
      dz2 = 1._RP/D%dz**2

      dkx_h = pi/(D%gnx)
      dky_h = pi/(D%gny)
      dkz_h = pi/(D%gnz)

      !$omp parallel private(i,j,k)

      !$omp workshare
      D%rwork = RHS
      !$omp end workshare


      !$omp end parallel
#ifdef MPI
      call PoisFFT_Plan_Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%rwork)
#endif
      !$omp parallel private(i,j,k)


      !$omp sections
      !$omp section
      forall(i=1:D%nx)  denomx(i) = (cos((i+D%offx)*dkx_h)-1.0_RP) * dx2 * 2
      !$omp section
      forall(j=1:D%ny)  denomy(j) = (cos((j+D%offy)*dky_h)-1.0_RP) * dy2 * 2
      !$omp section
      forall(k=1:D%nz)  denomz(k) = (cos((k+D%offz)*dkz_h)-1.0_RP) * dz2 * 2
      !$omp end sections

      !$omp do
      do k=1,D%nz
        do j=1,D%ny
          do i=1,D%nx
            D%rwork(i,j,k) = D%rwork(i,j,k) / (denomx(i) + denomy(j) + denomz(k))
          end do
        end do
      end do
      !$omp end do

      !$omp end parallel
#ifdef MPI
      call PoisFFT_Plan_Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp workshare
      Phi = D%rwork/(8*D%nx*D%ny*D%nz)
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullDirichletStag




    subroutine PoisFFT_Solver3D_FullNeumannStag(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      real(RP) :: dx2, dy2, dz2
      real(RP) :: dkx_h, dky_h, dkz_h !half of ussual dkx
      real(RP) :: denomx(D%nx), denomy(D%ny), denomz(D%nz)
      integer i,j,k
      real t1,t2

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2
      dz2 = 1._RP/D%dz**2

      dkx_h = pi/(D%gnx)
      dky_h = pi/(D%gny)
      dkz_h = pi/(D%gnz)

      !$omp parallel private(i,j,k)

      !$omp workshare
      D%rwork = RHS
      !$omp end workshare


      !$omp end parallel
#ifdef MPI
      call PoisFFT_Plan_Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp single
      if (D%offx==0.and.D%offy==0.and.D%offz==0) D%rwork(1,1,1) = 0
      !$omp end single

      !$omp sections
      !$omp section
      forall(i=1:D%nx)  denomx(i) = (cos((i-1+D%offx)*dkx_h)-1.0_RP) * dx2 * 2
      !$omp section
      forall(j=1:D%ny)  denomy(j) = (cos((j-1+D%offy)*dky_h)-1.0_RP) * dy2 * 2
      !$omp section
      forall(k=1:D%nz)  denomz(k) = (cos((k-1+D%offz)*dkz_h)-1.0_RP) * dz2 * 2
      !$omp end sections

      if (D%offx==0.and.D%offy==0.and.D%offz==0) then
#define xwork rwork
#include "loop_nest_3d.f90"
#undef xwork
        !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
        ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
      
        !$omp single
        D%rwork(1,1,1) = 0
        !$omp end single
      else
        !$omp do
        do k=1,D%nz
          do j=1,D%ny
            do i=1,D%nx
              D%rwork(i,j,k) = D%rwork(i,j,k) / (denomx(i) + denomy(j) + denomz(k))
            end do
          end do
        end do
        !$omp end do
      end if


      !$omp end parallel
#ifdef MPI
      call PoisFFT_Plan_Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp workshare
      Phi = D%rwork/(8*D%gnx*D%gny*D%gnz)
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullNeumannStag




    subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      real(RP) :: dx2, dy2, dz2
      real(RP) :: dkx, dky, dkz
      real(RP) :: denomx(D%nx), denomy(D%ny), denomz(D%nz)
      integer i,j,k
      integer tid      !thread id

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2
      dz2 = 1._RP/D%dz**2

      dkx = 2._RP*pi/(D%nx)
      dky = 2._RP*pi/(D%ny)
      dkz = 2._RP*pi/(D%nz)

      !$omp parallel private(tid,i,j,k)
      tid = 1
      !$ tid = omp_get_thread_num()+1

      ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
      !$omp do
      do j=1,D%ny
        do i=1,D%nx
          D%Solvers1D(tid)%rwork = RHS(i,j,:)


          call Execute(D%Solvers1D(tid)%forward,D%Solvers1D(tid)%rwork)


          Phi(i,j,:) = D%Solvers1D(tid)%rwork
        end do
      end do
      !$omp end do

      !$omp sections
      !$omp section
      forall(i=1:D%nx)  denomx(i) = (cos((i-1)*dkx)-1.0_RP) * dx2 * 2
      !$omp section
      forall(j=1:D%ny)  denomy(j) = (cos((j-1)*dky)-1.0_RP) * dy2 * 2
      !$omp section
      forall(k=1:D%nz)  denomz(k) = (cos((k-1)*pi/D%nz)-1.0_RP) * dz2 * 2
      !$omp end sections

      !$omp do
      do k=1,D%nz
        D%Solvers2D(tid)%cwork(1:D%nx,1:D%ny) = cmplx(Phi(1:D%nx,1:D%ny,k),0._RP,CP)


        call Execute(D%Solvers2D(tid)%forward, D%Solvers2D(tid)%cwork)


        if (k==1) then
          do j=2,D%ny
            do i=2,D%nx
              D%Solvers2D(tid)%cwork(i,j) = D%Solvers2D(tid)%cwork(i,j)&
                                              / (denomx(i) + denomy(j))
            end do
          end do
          do i=2,D%nx
              D%Solvers2D(tid)%cwork(i,1) = D%Solvers2D(tid)%cwork(i,1)&
                                              / (denomx(i))
          end do
          do j=2,D%ny
              D%Solvers2D(tid)%cwork(1,j) = D%Solvers2D(tid)%cwork(1,j)&
                                              / (denomy(j))
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          D%Solvers2D(tid)%cwork(1,1) = 0

        else

          do j=1,D%ny
            do i=1,D%nx
              D%Solvers2D(tid)%cwork(i,j) = D%Solvers2D(tid)%cwork(i,j)&
                                              / (denomx(i) + denomy(j) + denomz(k))
            end do
          end do

        endif


        call Execute(D%Solvers2D(tid)%backward, D%Solvers2D(tid)%cwork)


        Phi(:,:,k) = real(D%Solvers2D(tid)%cwork,RP) / (2 * D%nx * D%ny * D%nz)


      end do
      !$omp end do


      !$omp do
      do j=1,D%ny
        do i=1,D%nx
          D%Solvers1D(tid)%rwork = Phi(i,j,:)


          call Execute(D%Solvers1D(tid)%backward,D%Solvers1D(tid)%rwork)


          Phi(i,j,:) = D%Solvers1D(tid)%rwork
        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine PoisFFT_Solver3D_PPNs





    subroutine PoisFFT_Solver3D_Init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer i

      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel


      if (all(D%BCs==PoisFFT_Periodic)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call data_allocate_complex(D)

        D%forward = PoisFFT_Plan3D_QuickCreate(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan3D_QuickCreate(D, [FFT_Complex, FFTW_BACKWARD])

      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan3D_QuickCreate(D, [(FFT_RealOdd10, i=1,3)])
        D%backward = PoisFFT_Plan3D_QuickCreate(D, [(FFT_RealOdd01, i=1,3)])

      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan3D_QuickCreate(D, [(FFT_RealEven10, i=1,3)])
        D%backward = PoisFFT_Plan3D_QuickCreate(D, [(FFT_RealEven01, i=1,3)])

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. all(D%BCs(5:6)==PoisFFT_NeumannStag)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(1)

        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call data_allocate_real(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D_QuickCreate(D%Solvers1D(1), [FFT_RealEven10])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D_QuickCreate(D%Solvers1D(1), [FFT_RealEven01])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          call data_allocate_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call data_allocate_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D_QuickCreate(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D_QuickCreate(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          call data_allocate_complex(D%Solvers2D(i))
        end do

      endif
    end subroutine PoisFFT_Solver3D_Init























    subroutine PoisFFT_Solver2D_New(D,nx,ny,dx,dy,BCs)
      type(PoisFFT_Solver2D), intent(out) :: D

      integer, intent(in)   :: nx, ny
      real(RP), intent(in)  :: dx, dy
      integer, intent(in)   :: bcs(4)

      D%dx = dx
      D%dy = dy

      D%nx = nx
      D%ny = ny

      D%gnx = nx
      D%gny = ny

      D%cnt = D%nx * D%ny

      D%BCs = BCs

      !create fftw plans and allocate working array
      call Init(D)
    end subroutine PoisFFT_Solver2D_New



    subroutine PoisFFT_Solver2D_Execute(D,Phi,RHS)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:)
      real(RP), intent(in)  :: RHS(:,:)

      integer   :: ngPhi(2), ngRHS(2)


      ngPhi = (ubound(Phi)-[D%nx,D%ny])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny])/2


!       if (all(D%BCs==PoisFFT_Periodic)) then
!
!         call PoisFFT_Solver2D_FullPeriodic(D,&
!                  Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
!                      ngPhi(2)+1:ngPhi(2)+D%ny),&
!                  RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
!                      ngRHS(2)+1:ngRHS(2)+D%ny))
!
!       else if (all(D%BCs==PoisFFT_DirichletStag)) then
!
!         call PoisFFT_Solver2D_FullDirichletStag(D,&
!                  Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
!                      ngPhi(2)+1:ngPhi(2)+D%ny),&
!                  RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
!                      ngRHS(2)+1:ngRHS(2)+D%ny))

      if (all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver2D_FullNeumannStag(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny))
      endif

    end subroutine PoisFFT_Solver2D_Execute





    subroutine PoisFFT_Solver2D_FullNeumannStag(D, Phi, RHS)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:)
      real(RP), intent(in)  :: RHS(:,:)
      real(RP) :: dx2, dy2
      real(RP) :: dkx_h, dky_h
      real(RP) :: denomx(D%nx), denomy(D%ny)
      integer i,j

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2

      dkx_h = pi/(D%nx)
      dky_h = pi/(D%ny)

      ! Forward FFT of RHS
      D%rwork = RHS


      call Execute(D%forward, D%rwork)


      forall(i=1:D%nx)  denomx(i) = (cos((i-1)*dkx_h)-1.0_RP) * dx2 * 2

      forall(j=1:D%ny)  denomy(j) = (cos((j-1)*dky_h)-1.0_RP) * dy2 * 2

      D%rwork(1,1) = 0

      do j = 1,D%ny
        do i = 1,D%nx
          if (i==1.and.j==1) cycle
          D%rwork(i,j) = D%rwork(i,j) / (denomx(i) + denomy(j))
        end do
      end do

      call Execute(D%backward, D%rwork)

      Phi = D%rwork/(4*D%nx*D%ny)

    end subroutine PoisFFT_Solver2D_FullNeumannStag



    subroutine Poisfft_Solver2D_Init(D)
      type(PoisFFT_Solver2D), intent(inout) :: D
      integer i


      if (all(D%BCs==PoisFFT_Periodic)) then

        D%forward = PoisFFT_Plan2D_QuickCreate(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan2D_QuickCreate(D, [FFT_Complex, FFTW_BACKWARD])

        call data_allocate_complex(D)


      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        D%forward = PoisFFT_Plan2D_QuickCreate(D, [(FFT_RealOdd10, i=1,2)])
        D%backward = PoisFFT_Plan2D_QuickCreate(D, [(FFT_RealOdd01, i=1,2)])

        call data_allocate_real(D)


      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        D%forward = PoisFFT_Plan2D_QuickCreate(D, [(FFT_RealEven10, i=1,2)])
        D%backward = PoisFFT_Plan2D_QuickCreate(D, [(FFT_RealEven01, i=1,2)])

        call data_allocate_real(D)


      endif
    end subroutine PoisFFT_Solver2D_Init






















    subroutine PoisFFT_Solver1D_New(D,nx,dx,BCs)
      type(PoisFFT_Solver1D), intent(out) :: D

      integer, intent(in)   :: nx
      real(RP), intent(in)  :: dx
      integer, intent(in)   :: bcs(2)

      D%dx = dx

      D%nx = nx

      D%gnx = nx

      D%cnt = D%nx

      D%BCs = BCs

      !create fftw plans and allocate working array
      call Init(D)
    end subroutine PoisFFT_Solver1D_New



    subroutine PoisFFT_Solver1D_Execute(D,Phi,RHS)
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

      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call PoisFFT_Solver1D_FullDirichletStag(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))

      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver1D_FullNeumannStag(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx))
      endif

    end subroutine PoisFFT_Solver1D_Execute





    subroutine PoisFFT_Solver1D_FullPeriodic(D, Phi, RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      real(RP) :: dx2
      real(RP) :: dkx
      integer i

      dx2 = 1._RP/D%dx**2

      dkx = 2._RP*pi/(D%nx)

      ! Forward FFT of RHS
      D%cwork = cmplx(RHS,0._RP,CP)


      call Execute(D%forward, D%cwork)


      D%cwork(1) = 0
      forall(i=2:D%nx) &
        D%cwork(i) = D%cwork(i) / (2.0_RP * (cos((i-1)*dkx)-1.0_RP)*dx2)

      call Execute(D%backward, D%cwork)

      Phi = real(D%cwork,RP) / (D%nx)

    end subroutine PoisFFT_Solver1D_FullPeriodic





    subroutine PoisFFT_Solver1D_FullDirichletStag(D, Phi, RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      real(RP) :: dx2
      real(RP) :: dkx_h !half of ussual dkx
      integer i

      dx2 = 1._RP/D%dx**2
      dkx_h = pi/(D%nx)

      ! Forward FFT of RHS
      D%rwork = RHS


      call Execute(D%forward, D%rwork)


      forall(i=1:D%nx) &
        D%rwork(i) = D%rwork(i) / (2.0_RP * (cos((i)*dkx_h)-1.0_RP)*dx2)

      call Execute(D%backward, D%rwork)

      Phi = D%rwork/(2*D%nx)

    end subroutine PoisFFT_Solver1D_FullDirichletStag




    subroutine PoisFFT_Solver1D_FullNeumannStag(D, Phi, RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      real(RP) :: dx2
      real(RP) :: dkx_h
      integer i

      dx2 = 1._RP/D%dx**2

      dkx_h = pi/(D%nx)

      ! Forward FFT of RHS
      D%rwork = RHS


      call Execute(D%forward, D%rwork)


      D%rwork = D%rwork/(D%nx)

      D%rwork(1) = 0
      forall(i=2:D%nx) &
        D%rwork(i) = D%rwork(i) / (2.0_RP * (cos((i-1)*dkx_h)-1.0_RP)*dx2)

      call Execute(D%backward, D%rwork)

      Phi = D%rwork/2

    end subroutine PoisFFT_Solver1D_FullNeumannStag



    subroutine Poisfft_Solver1D_Init(D)
      type(PoisFFT_Solver1D), intent(inout) :: D
      integer i


      if (all(D%BCs==PoisFFT_Periodic)) then

        D%forward = PoisFFT_Plan1D_QuickCreate(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan1D_QuickCreate(D, [FFT_Complex, FFTW_BACKWARD])

        call data_allocate_complex(D)


      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        D%forward = PoisFFT_Plan1D_QuickCreate(D, [FFT_RealOdd10])
        D%backward = PoisFFT_Plan1D_QuickCreate(D, [FFT_RealOdd01])

        call data_allocate_real(D)


      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        D%forward = PoisFFT_Plan1D_QuickCreate(D, [FFT_RealEven10])
        D%backward = PoisFFT_Plan1D_QuickCreate(D, [FFT_RealEven01])

        call data_allocate_real(D)


      endif
    end subroutine PoisFFT_Solver1D_Init
