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
    module procedure PoisFFT_Solver1D_New
  end interface

  interface PoisFFT_Solver2D
    module procedure PoisFFT_Solver2D_New
  end interface

  interface PoisFFT_Solver3D
    module procedure PoisFFT_Solver3D_New
  end interface

  interface Finalize
    module procedure PoisFFT_Solver1D_Finalize
    module procedure PoisFFT_Solver2D_Finalize
    module procedure PoisFFT_Solver3D_Finalize
  end interface Finalize


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



    function PoisFFT_Solver3D_New(nxyz,dxyz,BCs,gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver3D) :: D

      integer, intent(in)   :: nxyz(3)
      real(RP), intent(in)  :: dxyz(3)
      integer, intent(in)   :: bcs(6)
      integer, intent(in), optional :: gnxyz(3)
      integer, intent(in), optional :: offs(3)
      integer, intent(in), optional  :: mpi_comm
      integer, intent(in), optional :: nthreads

      D%nxyz = nxyz

      D%dx = dxyz(1)
      D%dy = dxyz(2)
      D%dz = dxyz(3)

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
      
#ifdef MPI      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
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
    end function PoisFFT_Solver3D_New


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
      
#endif

      if (direction==1) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%BCs = D3D%BCs(1:2)
      else if (direction==2) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%BCs = D3D%BCs(3:4)
      else
        D%dx = D3D%dz
        D%nx = D3D%nz
        D%gnx = D3D%gnz
        D%offx = D3D%offz
        D%BCs = D3D%BCs(5:6)
      endif
      
      D%nxyz = [D%nx]

      D%cnt = D%nx
      
      D%gcnt = int(D%gnx, kind(D%gcnt))
    end function PoisFFT_Solver1D_From3D


    function PoisFFT_Solver1D_Many_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver1D_Many)        :: D
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
      
#endif

      if (direction==1) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%BCs = D3D%BCs(1:2)
        D%howmany = D3D%ny * D3D%nz
        D%workdims = [D3D%nx, D3D%ny, D3D%nz]
      else if (direction==2) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%BCs = D3D%BCs(3:4)
        D%howmany = D3D%nx * D3D%nz
        D%workdims = [D3D%ny, D3D%nx, D3D%nz]
      else
        D%dx = D3D%dz
        D%nx = D3D%nz
        D%gnx = D3D%gnz
        D%offx = D3D%offz
        D%BCs = D3D%BCs(5:6)
        D%howmany = D3D%nx * D3D%ny
        D%workdims = [D3D%nz, D3D%ny, D3D%nx]
      endif

      D%nxyz = [D%nx]
        
      D%cnt = D%nx
      
      D%gcnt = int(D%gnx, kind(D%gcnt))
    end function PoisFFT_Solver1D_Many_From3D


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
        else if (direction==2) then
          call MPI_Cart_sub(D3D%mpi%comm, [.true.,.false.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
        else
          call MPI_Cart_sub(D3D%mpi%comm, [.false.,.true.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
        end if
      else
        stop "Not implemented."
      end if
      
#endif

      if (direction==1) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%offy = D3D%offz
        D%BCs = D3D%BCs(3:6)
      else if (direction==2) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%offy = D3D%offz
        D%BCs = D3D%BCs([1,2,5,6])
      else
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%dy = D3D%dy
        D%ny = D3D%ny
        D%gny = D3D%gny
        D%offy = D3D%offy
        D%BCs = D3D%BCs(1:4)
      endif

      D%nxyz = [D%nx, D%ny]

      D%cnt = D%nx * D%ny
      
      D%gcnt = int(D%gnx, kind(D%gcnt)) * int(D%gny, kind(D%gcnt))
    end function PoisFFT_Solver2D_From3D


    function PoisFFT_Solver2D_Many_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver2D_Many)        :: D
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
        else if (direction==2) then
          call MPI_Cart_sub(D3D%mpi%comm, [.true.,.false.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
        else
          call MPI_Cart_sub(D3D%mpi%comm, [.false.,.true.], D%mpi%comm, ie)
          if (ie/=0) stop "Error executing MPI_Cart_sub."
        end if
      else
        stop "Not implemented."
      end if
      
#endif

      if (direction==1) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%gnx = D3D%gny
        D%offx = D3D%offy
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%offy = D3D%offz
        D%BCs = D3D%BCs(3:6)
        D%howmany = D3D%nx
        D%workdims = [D3D%ny, D3D%nz, D3D%nx]
      else if (direction==2) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%gny = D3D%gnz
        D%offy = D3D%offz
        D%BCs = D3D%BCs([1,2,5,6])
        D%howmany = D3D%ny
        D%workdims = [D3D%nx, D3D%nz, D3D%ny]
      else
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%gnx = D3D%gnx
        D%offx = D3D%offx
        D%dy = D3D%dy
        D%ny = D3D%ny
        D%gny = D3D%gny
        D%offy = D3D%offy
        D%BCs = D3D%BCs(1:4)
        D%howmany = D3D%nz
        D%workdims = [D3D%nx, D3D%ny, D3D%nz]
      endif
      
      D%nxyz = [D%nx, D%ny]

      D%cnt = D%nx * D%ny
      
      D%gcnt = int(D%gnx, kind(D%gcnt)) * int(D%gny, kind(D%gcnt))
    end function PoisFFT_Solver2D_Many_From3D


    subroutine PoisFFT_Solver1D_Finalize(D)
      type(PoisFFT_Solver1D), intent(inout) :: D

      call Destroy(D%forward)
      call Destroy(D%backward)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

    endsubroutine PoisFFT_Solver1D_Finalize



    subroutine PoisFFT_Solver2D_Finalize(D)
      type(PoisFFT_Solver2D), intent(inout) :: D

      call Destroy(D%forward)
      call Destroy(D%backward)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

    endsubroutine PoisFFT_Solver2D_Finalize



    subroutine PoisFFT_Solver3D_Finalize(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: i

      call Destroy(D%forward)
      call Destroy(D%backward)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

      if (allocated(D%Solvers1D)) then
        do i=1,size(D%Solvers1D)
          call Finalize(D%Solvers1D(i))
        end do
      endif

      if (allocated(D%Solvers2D)) then
        do i=1,size(D%Solvers2D)
          call Finalize(D%Solvers2D(i))
        end do
      endif

!       if (allocated(D%Solver1D)) call Finalize(D%Solver1D)
! 
!       if (allocated(D%Solver2D)) call Finalize(D%Solver2D)

    endsubroutine PoisFFT_Solver3D_Finalize



    subroutine PoisFFT_Solver3D_FullPeriodic(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k

      !$omp parallel private(i,j,k)

      !$omp workshare
      D%cwork = cmplx(RHS,0._RP,CP)
      !$omp end workshare

      !$omp end parallel
#ifdef MPI
      call Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%cwork)
#endif
      !$omp parallel private(i,j,k)

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
              D%cwork(i,j,k) = D%cwork(i,j,k) / (D%denomx(i) + D%denomy(j) + D%denomz(k))
            end do
          end do
        end do
        !$omp end do
      end if

      !$omp end parallel
#ifdef MPI
      call Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%cwork)
#endif
      !$omp parallel

      !$omp workshare
       Phi = real(D%cwork,RP) / (D%gcnt)
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullPeriodic




    subroutine PoisFFT_Solver3D_FullDirichletStag(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k

      !$omp parallel private(i,j,k)

      !$omp workshare
      D%rwork = RHS
      !$omp end workshare


      !$omp end parallel
#ifdef MPI
      call Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp do
      do k=1,D%nz
        do j=1,D%ny
          do i=1,D%nx
            D%rwork(i,j,k) = D%rwork(i,j,k) / (D%denomx(i) + D%denomy(j) + D%denomz(k))
          end do
        end do
      end do
      !$omp end do

      !$omp end parallel
#ifdef MPI
      call Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp workshare
      Phi = D%rwork/(8*D%gcnt)
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullDirichletStag




    subroutine PoisFFT_Solver3D_FullNeumannStag(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k

      !$omp parallel private(i,j,k)

      !$omp workshare
      D%rwork = RHS
      !$omp end workshare


      !$omp end parallel
#ifdef MPI
      call Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp single
      if (D%offx==0.and.D%offy==0.and.D%offz==0) D%rwork(1,1,1) = 0
      !$omp end single

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
              D%rwork(i,j,k) = D%rwork(i,j,k) / (D%denomx(i) + D%denomy(j) + D%denomz(k))
            end do
          end do
        end do
        !$omp end do
      end if


      !$omp end parallel
#ifdef MPI
      call Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%rwork)
#endif
      !$omp parallel private(i,j,k)

      !$omp workshare
      Phi = D%rwork/(8*D%gcnt)
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullNeumannStag



#ifdef MPI

#ifdef MANY_VARIANT
!     subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
!       type(PoisFFT_Solver3D), intent(inout) :: D
!       real(RP), intent(out) :: Phi(:,:,:)
!       real(RP), intent(in)  :: RHS(:,:,:)
!       integer i,j,k
! 
!       interface
!         subroutine MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE, &
!                                  RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE, COMM, IERROR)
!           import
!           REAL(RP)    SENDBUF(*), RECVBUF(*)
!           INTEGER    SENDCOUNTS(*), SDISPLS(*), SENDTYPE
!           INTEGER    RECVCOUNTS(*), RDISPLS(*), RECVTYPE
!           INTEGER    COMM, IERROR
!         end subroutine
!       end interface
!       integer :: l, ie
! 
! 
! 
! 
!       if (D%Solver1D%mpi_transpose_needed) then
!       
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               D%mpi%tmp1(k,j,i) = RHS(i,j,k)
!             end do
!           end do
!         end do
! 
!         
!         !step2 exchange
!         call MPI_AllToAllV(D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
!                            D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
!                            D%Solver1D%mpi%comm, ie)
!         !step3 local reordering of blocks
! 
!         do l = 1, D%mpi%np
!           do k = 0, D%mpi%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, D%mpi%rnzs(l)-1
!                 D%Solver1D%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1) = &
!                   D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l))
!               end do
!             end do
!           end do
!         end do
!      
!         ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
!         call Execute(D%Solver1D%forward,D%Solver1D%rwork)
! 
!         
! 
!         !step3' local reordering of blocks
!         do l = 1, D%mpi%np
!           do k = 0, D%mpi%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, D%mpi%rnzs(l)-1
!                 D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l)) = &
!                   D%Solver1D%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1)
!               end do
!             end do
!           end do
!         end do
! 
!       
!         !step2' exchange
!         call MPI_AllToAllV(D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
!                            D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
!                            D%Solver1D%mpi%comm, ie)
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               D%Solver2D%cwork(i,j,k) = D%mpi%tmp1(k,j,i)
!             end do
!           end do
!         end do
!       else
!         ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
!         do j = 1, D%ny
!           do k = 1, D%nz
!             do i = 1, D%nx
!               D%Solver1D%rwork(k,j,i) = RHS(i,j,k)
!             end do
!           end do
!         end do
! 
!         call Execute(D%Solver1D%forward, D%Solver1D%rwork)
! 
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               D%Solver2D%cwork(i,j,k) = D%Solver1D%rwork(k,j,i)
!             end do
!           end do
!         end do
!       end if
! 
!       call Execute_MPI(D%Solver2D%forward)
! 
!       do k=1,D%nz
!         if (k==1.and.D%offz==0) then
!           do j=2,D%ny
!             do i=2,D%nx
!               D%Solver2D%cwork(i,j,k) = D%Solver2D%cwork(i,j,k)&
!                                               / (D%denomx(i) + D%denomy(j))
!             end do
!           end do
!           do i=2,D%nx
!               D%Solver2D%cwork(i,1,k) = D%Solver2D%cwork(i,1,k)&
!                                               / (D%denomx(i))
!           end do
!           do j=2,D%ny
!               D%Solver2D%cwork(1,j,k) = D%Solver2D%cwork(1,j,k)&
!                                               / (D%denomy(j))
!           end do
!           !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
!           ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
!           D%Solver2D%cwork(1,1,1) = 0
! 
!         else
! 
!           do j=1,D%ny
!             do i=1,D%nx
!               D%Solver2D%cwork(i,j,k) = D%Solver2D%cwork(i,j,k)&
!                                               / (D%denomx(i) + D%denomy(j) + D%denomz(k))
!             end do
!           end do
! 
!         endif
!       end do
! 
!       call Execute_MPI(D%Solver2D%backward)
! 
! 
!       if (D%Solver1D%mpi_transpose_needed) then
! 
!         !step1 local transpose
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               D%mpi%tmp1(k,j,i) = real(D%Solver2D%cwork(i,j,k),RP) / (2 * D%gcnt)
!             end do
!           end do
!         end do
! 
!         !step2 exchange
!         call MPI_AllToAllV(D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
!                            D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
!                            D%Solver1D%mpi%comm, ie)
! 
!         !step3 local reordering of blocks
!         do l = 1, D%mpi%np
!           do k = 0, D%mpi%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, D%mpi%rnzs(l)-1
!                 D%Solver1D%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1) = &
!                   D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l))
!               end do
!             end do
!           end do
!         end do
!       
!         ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
!         call Execute(D%Solver1D%backward, D%Solver1D%rwork)
! 
!         !step3' local reordering of blocks
!         do l = 1, D%mpi%np
!           do k = 0, D%mpi%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, D%mpi%rnzs(l)-1
!                 D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l)) = &
!                   D%Solver1D%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1)
!               end do
!             end do
!           end do
!         end do
!       
!         !step2' exchange
!         call MPI_AllToAllV(D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
!                            D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
!                            D%Solver1D%mpi%comm, ie)
! 
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               Phi(i,j,k) = D%mpi%tmp1(k,j,i)
!             end do
!           end do
!         end do
!       else
!         ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
! 
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               D%Solver1D%rwork(k,j,i) = real(D%Solver2D%cwork(i,j,k),RP) / (2 * D%gcnt)
!             end do
!           end do
!         end do
! 
!         call Execute(D%Solver1D%backward, D%Solver1D%rwork)
! 
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               Phi(i,j,k) = D%Solver1D%rwork(k,j,i)
!             end do
!           end do
!         end do
! 
!       end if
!     end subroutine PoisFFT_Solver3D_PPNs

!MANY_VARIANT
#else

    subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k

      interface
        subroutine MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE, &
                                 RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE, COMM, IERROR)
          import
          REAL(RP)    SENDBUF(*), RECVBUF(*)
          INTEGER    SENDCOUNTS(*), SDISPLS(*), SENDTYPE
          INTEGER    RECVCOUNTS(*), RDISPLS(*), RECVTYPE
          INTEGER    COMM, IERROR
        end subroutine
      end interface
      integer :: l, ie

      if (D%Solvers1D(1)%mpi_transpose_needed) then
      
        do k = 1, D%nz
          do j = 1, D%ny
            do i = 1, D%nx
              D%mpi%tmp1(k,j,i) = RHS(i,j,k)
            end do
          end do
        end do

        
        !step2 exchange
        call MPI_AllToAllV(D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
                           D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
                           D%Solvers1D(1)%mpi%comm, ie)
        !step3 local reordering of blocks

        do l = 1, D%mpi%np
          do k = 0, D%mpi%rnxs(1)-1
            do j = 0, D%ny-1
              do i = 0, D%mpi%rnzs(l)-1
                D%mpi%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1) = &
                  D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l))
              end do
            end do
          end do
        end do
     
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        do k=1,size(D%mpi%rwork,3)
          do j=1,size(D%mpi%rwork,2)
            call Execute(D%Solvers1D(1)%forward,D%mpi%rwork(:,j,k))
          end do
        end do
        

        !step3' local reordering of blocks
        do l = 1, D%mpi%np
          do k = 0, D%mpi%rnxs(1)-1
            do j = 0, D%ny-1
              do i = 0, D%mpi%rnzs(l)-1
                D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l)) = &
                  D%mpi%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1)
              end do
            end do
          end do
        end do

      
        !step2' exchange
        call MPI_AllToAllV(D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
                           D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
                           D%Solvers1D(1)%mpi%comm, ie)
        do k = 1, D%nz
          do j = 1, D%ny
            do i = 1, D%nx
              Phi(i,j,k) = D%mpi%tmp1(k,j,i)
            end do
          end do
        end do
      else
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(1)%rwork = RHS(i,j,:)


            call Execute(D%Solvers1D(1)%forward,D%Solvers1D(1)%rwork)


            Phi(i,j,:) = D%Solvers1D(1)%rwork
          end do
        end do
      end if
      

      do k=1,D%nz
        D%Solvers2D(1)%cwork(1:D%nx,1:D%ny) = cmplx(Phi(1:D%nx,1:D%ny,k),0._RP,CP)

        call Execute_MPI(D%Solvers2D(1)%forward)

        if (k==1.and.D%offz==0) then
          do j = 1,D%ny
            do i = max(3-j-D%offx-D%offy,1),D%nx
              D%Solvers2D(1)%cwork(i,j) = D%Solvers2D(1)%cwork(i,j) / (D%denomx(i) + D%denomy(j))
            end do
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          if (D%offx==0.and.D%offy==0) D%Solvers2D(1)%cwork(1,1) = 0

        else

          do j=1,D%ny
            do i=1,D%nx
              D%Solvers2D(1)%cwork(i,j) = D%Solvers2D(1)%cwork(i,j)&
                                              / (D%denomx(i) + D%denomy(j) + D%denomz(k))
            end do
          end do

        endif

        call Execute_MPI(D%Solvers2D(1)%backward)

        Phi(:,:,k) = real(D%Solvers2D(1)%cwork,RP) / (2 * D%gcnt)


      end do


      if (D%Solvers1D(1)%mpi_transpose_needed) then

        !step1 local transpose
        do k = 1, D%nz
          do j = 1, D%ny
            do i = 1, D%nx
              D%mpi%tmp1(k,j,i) = Phi(i,j,k)
            end do
          end do
        end do

        !step2 exchange
        call MPI_AllToAllV(D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
                           D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
                           D%Solvers1D(1)%mpi%comm, ie)

        !step3 local reordering of blocks
        do l = 1, D%mpi%np
          do k = 0, D%mpi%rnxs(1)-1
            do j = 0, D%ny-1
              do i = 0, D%mpi%rnzs(l)-1
                D%mpi%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1) = &
                  D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l))
              end do
            end do
          end do
        end do
      
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        do k=1,size(D%mpi%rwork,3)
          do j=1,size(D%mpi%rwork,2)
            call Execute(D%Solvers1D(1)%backward, D%mpi%rwork(:,j,k))
          end do
        end do

        !step3' local reordering of blocks
        do l = 1, D%mpi%np
          do k = 0, D%mpi%rnxs(1)-1
            do j = 0, D%ny-1
              do i = 0, D%mpi%rnzs(l)-1
                D%mpi%tmp2(i + j*D%mpi%rnzs(l) + k*(D%ny*D%mpi%rnzs(l)) + D%mpi%rdispls(l)) = &
                  D%mpi%rwork(i+sum(D%mpi%rnzs(1:l-1))+1, j+1, k+1)
              end do
            end do
          end do
        end do
      
        !step2' exchange
        call MPI_AllToAllV(D%mpi%tmp2, D%mpi%rcounts, D%mpi%rdispls, MPI_RP, &
                           D%mpi%tmp1, D%mpi%scounts, D%mpi%sdispls, MPI_RP, &
                           D%Solvers1D(1)%mpi%comm, ie)

        do k = 1, D%nz
          do j = 1, D%ny
            do i = 1, D%nx
              Phi(i,j,k) = D%mpi%tmp1(k,j,i)
            end do
          end do
        end do
      else
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(1)%rwork = Phi(i,j,:)


            call Execute(D%Solvers1D(1)%backward,D%Solvers1D(1)%rwork)


            Phi(i,j,:) = D%Solvers1D(1)%rwork
          end do
        end do
      end if
    end subroutine PoisFFT_Solver3D_PPNs
#endif
    
#elif defined(MANY_VARIANT)

!     subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
!       type(PoisFFT_Solver3D), intent(inout) :: D
!       real(RP), intent(out) :: Phi(:,:,:)
!       real(RP), intent(in)  :: RHS(:,:,:)
!       integer i,j,k
!       integer tid      !thread id
! 
! 
!       !$omp parallel private(tid,i,j,k)
!       tid = 1
!       !$ tid = omp_get_thread_num()+1
! 
!       ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
!       !$omp do
!       do j = 1, D%ny
!         do k = 1, D%nz
!           do i = 1, D%nx
!             D%Solver1D%rwork(k,j,i) = RHS(i,j,k)
!           end do
!         end do
!       end do
!       !$omp end do
!       !$omp single
!       call Execute(D%Solver1D%forward, D%Solver1D%rwork)
!       !$omp end single
!       !$omp do
!       do k = 1, D%nz
!         do j = 1, D%ny
!           do i = 1, D%nx
!             D%Solver2D%cwork(i,j,k) = D%Solver1D%rwork(k,j,i)
!           end do
!         end do
!       end do
!       !$omp end do
! 
! 
!       !$omp single
!       call Execute(D%Solver2D%forward, D%Solver2D%cwork)
!       !$omp end single
!       
!       !$omp do
!       do k=1,D%nz
!         if (k==1.and.D%offz==0) then
!           do j=2,D%ny
!             do i=2,D%nx
!               D%Solver2D%cwork(i,j,k) = D%Solver2D%cwork(i,j,k)&
!                                               / (D%denomx(i) + D%denomy(j))
!             end do
!           end do
!           do i=2,D%nx
!               D%Solver2D%cwork(i,1,k) = D%Solver2D%cwork(i,1,k)&
!                                               / (D%denomx(i))
!           end do
!           do j=2,D%ny
!               D%Solver2D%cwork(1,j,k) = D%Solver2D%cwork(1,j,k)&
!                                               / (D%denomy(j))
!           end do
!           !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
!           ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
!           D%Solver2D%cwork(1,1,1) = 0
! 
!         else
! 
!           do j=1,D%ny
!             do i=1,D%nx
!               D%Solver2D%cwork(i,j,k) = D%Solver2D%cwork(i,j,k)&
!                                               / (D%denomx(i) + D%denomy(j) + D%denomz(k))
!             end do
!           end do
! 
!         endif
!       end do
!       !$omp end do
!       
!       !$omp single
!       call Execute(D%Solver2D%backward, D%Solver2D%cwork)
!       !$omp end single
! 
!       !$omp do
!       do k = 1, D%nz
!         do j = 1, D%ny
!           do i = 1, D%nx
!             D%Solver1D%rwork(k,j,i) = real(D%Solver2D%cwork(i,j,k),RP) / (2 * D%gcnt)
!           end do
!         end do
!       end do
!       !$omp end do
!       !$omp single
!       call Execute(D%Solver1D%backward, D%Solver1D%rwork)
!       !$omp end single
!       !$omp do
!       do k = 1, D%nz
!         do j = 1, D%ny
!           do i = 1, D%nx
!             Phi(i,j,k) = D%Solver1D%rwork(k,j,i)
!           end do
!         end do
!       end do
!       !$omp end do
!       !$omp end parallel
!     end subroutine PoisFFT_Solver3D_PPNs
    
#else

    subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k
      integer tid      !thread id

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
      
      
      !$omp do
      do k=1,D%nz
        D%Solvers2D(tid)%cwork(1:D%nx,1:D%ny) = cmplx(Phi(1:D%nx,1:D%ny,k),0._RP,CP)


        call Execute(D%Solvers2D(tid)%forward, D%Solvers2D(tid)%cwork)

        if (k==1) then
          do j=2,D%ny
            do i=2,D%nx
              D%Solvers2D(tid)%cwork(i,j) = D%Solvers2D(tid)%cwork(i,j)&
                                              / (D%denomx(i) + D%denomy(j))
            end do
          end do
          do i=2,D%nx
              D%Solvers2D(tid)%cwork(i,1) = D%Solvers2D(tid)%cwork(i,1)&
                                              / (D%denomx(i))
          end do
          do j=2,D%ny
              D%Solvers2D(tid)%cwork(1,j) = D%Solvers2D(tid)%cwork(1,j)&
                                              / (D%denomy(j))
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          D%Solvers2D(tid)%cwork(1,1) = 0

        else

          do j=1,D%ny
            do i=1,D%nx
              D%Solvers2D(tid)%cwork(i,j) = D%Solvers2D(tid)%cwork(i,j)&
                                              / (D%denomx(i) + D%denomy(j) + D%denomz(k))
            end do
          end do

        endif

        call Execute(D%Solvers2D(tid)%backward, D%Solvers2D(tid)%cwork)

        Phi(:,:,k) = real(D%Solvers2D(tid)%cwork,RP) / (2 * D%gcnt)


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
    
#endif




    subroutine PoisFFT_Solver3D_Init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: i, j, k
      real(RP) :: dkx, dky, dkz, dkx_h, dky_h, dkz_h
#ifdef MPI
      integer :: ie
      
      interface
        subroutine MPI_ALLTOALL(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,&
                                RECVTYPE, COMM, IERROR)
            INTEGER   SENDBUF(*), RECVBUF(*)
            INTEGER   SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
            INTEGER   COMM, IERROR
        end subroutine
      end interface
#endif

      !$omp parallel
      !$omp single
      !$ D%nthreads = omp_get_num_threads()
      !$omp end single
      !$omp end parallel

      dkx_h = pi/(D%gnx)
      dky_h = pi/(D%gny)
      dkz_h = pi/(D%gnz)
      
      dkx = 2 * dkx_h
      dky = 2 * dky_h
      dkz = 2 * dkz_h

      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))
      allocate(D%denomz(D%nz))
      
      if (all(D%BCs==PoisFFT_Periodic)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call data_allocate_complex(D)

        D%forward = PoisFFT_Plan3D_Create(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan3D_Create(D, [FFT_Complex, FFTW_BACKWARD])

        forall(i=1:D%nx)  D%denomx(i) = 2 * (cos((i-1+D%offx)*dkx)-1.0_RP) / D%dx**2
        forall(j=1:D%ny)  D%denomy(j) = 2 * (cos((j-1+D%offy)*dky)-1.0_RP) / D%dy**2
        forall(k=1:D%nz)  D%denomz(k) = 2 * (cos((k-1+D%offz)*dkz)-1.0_RP) / D%dz**2
        
      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan3D_Create(D, [(FFT_RealOdd10, i=1,3)])
        D%backward = PoisFFT_Plan3D_Create(D, [(FFT_RealOdd01, i=1,3)])

        forall(i=1:D%nx)  D%denomx(i) = 2 * (cos((i+D%offx)*dkx_h)-1.0_RP) / D%dx**2
        forall(j=1:D%ny)  D%denomy(j) = 2 * (cos((j+D%offy)*dky_h)-1.0_RP) / D%dy**2
        forall(k=1:D%nz)  D%denomz(k) = 2 * (cos((k+D%offz)*dkz_h)-1.0_RP) / D%dz**2
      
      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan3D_Create(D, [(FFT_RealEven10, i=1,3)])
        D%backward = PoisFFT_Plan3D_Create(D, [(FFT_RealEven01, i=1,3)])

        forall(i=1:D%nx)  D%denomx(i) = 2 * (cos((i-1+D%offx)*dkx_h)-1.0_RP) / D%dx**2
        forall(j=1:D%ny)  D%denomy(j) = 2 * (cos((j-1+D%offy)*dky_h)-1.0_RP) / D%dy**2
        forall(k=1:D%nz)  D%denomz(k) = 2 * (cos((k-1+D%offz)*dkz_h)-1.0_RP) / D%dz**2

      else if (all(D%BCs(1:4)==PoisFFT_Periodic) .and. all(D%BCs(5:6)==PoisFFT_NeumannStag)) then
      
#ifdef MANY_VARIANT

! #ifndef MPI
!         if (D%nthreads>1) call PoisFFT_InitThreads(D%nthreads)
! #endif
!         allocate(D%Solver1D)
!         allocate(D%Solver2D)
! 
!         D%Solver1D = PoisFFT_Solver1D_Many_From3D(D,3)
! 
!         D%Solver2D = PoisFFT_Solver2D_Many_From3D(D,3)
! 
! #ifdef MPI
!         call MPI_Comm_rank(D%Solver1D%mpi%comm, D%mpi%rank, ie)
!         call MPI_Comm_size(D%Solver1D%mpi%comm, D%mpi%np, ie)
!         allocate(D%mpi%snxs(D%mpi%np))
!         allocate(D%mpi%snzs(D%mpi%np))
!         allocate(D%mpi%rnxs(D%mpi%np))
!         allocate(D%mpi%rnzs(D%mpi%np))
!         allocate(D%mpi%sdispls(D%mpi%np))
!         allocate(D%mpi%scounts(D%mpi%np))
!         allocate(D%mpi%rdispls(D%mpi%np))
!         allocate(D%mpi%rcounts(D%mpi%np))
!         
!         D%mpi%snxs(1:D%mpi%np-1) = (D%nx / D%mpi%np)
!         D%mpi%snxs(D%mpi%np) = D%nx - sum(D%mpi%snxs(1:D%mpi%np-1))
!         D%mpi%snzs = D%nz
!         D%mpi%scounts = D%mpi%snxs*D%ny*D%nz
!         D%mpi%sdispls(1) = 0
!         do i = 2, D%mpi%np
!           D%mpi%sdispls(i) = D%mpi%sdispls(i-1) + D%mpi%scounts(i-1)
!         end do
! 
!         call MPI_AllToAll(D%mpi%snxs, 1, MPI_INTEGER, &
!                           D%mpi%rnxs, 1, MPI_INTEGER, &
!                           D%Solver1D%mpi%comm, ie)
!         if (.not.all(D%mpi%rnxs(2:D%mpi%np)==D%mpi%rnxs(1))) stop "????"
!         call MPI_AllToAll(D%mpi%snzs, 1, MPI_INTEGER, &
!                           D%mpi%rnzs, 1, MPI_INTEGER, &
!                           D%Solver1D%mpi%comm, ie)
!         call MPI_AllToAll(D%mpi%scounts, 1, MPI_INTEGER, &
!                           D%mpi%rcounts, 1, MPI_INTEGER, &
!                           D%Solver1D%mpi%comm, ie)
! 
!         call MPI_AllToAll(D%mpi%sdispls, 1, MPI_INTEGER, &
!                           D%mpi%rdispls, 1, MPI_INTEGER, &
!                           D%Solver1D%mpi%comm, ie)
! 
!         D%mpi%rdispls(1) = 0
!         do i = 2, D%mpi%np
!           D%mpi%rdispls(i) = D%mpi%rdispls(i-1) + D%mpi%rcounts(i-1)
!         end do
!         !step1 local transpose
!         allocate(D%mpi%tmp1(1:D%nz,1:D%ny,1:D%nx))
!         allocate(D%mpi%tmp2(0:sum(D%mpi%rcounts)-1))
!         
!         D%Solver1D%workdims = [sum(D%mpi%rnzs), D%ny, D%mpi%rnxs(1)]
!         D%Solver1D%howmany = D%ny * D%mpi%rnxs(1)
! #endif
! 
! 
! 
! 
!         call data_allocate_real(D%Solver1D)
! 
!         D%Solver1D%forward = PoisFFT_Plan1D_Many_Create(D%Solver1D, [FFT_RealEven10])
!         D%Solver1D%backward = PoisFFT_Plan1D_Many_Create(D%Solver1D, [FFT_RealEven01])
! 
!         call data_allocate_complex(D%Solver2D)
! 
!         D%Solver2D%forward = PoisFFT_Plan2D_Many_Create(D%Solver2D, [FFT_Complex, FFTW_FORWARD])
!         D%Solver2D%backward = PoisFFT_Plan2D_Many_Create(D%Solver2D, [FFT_Complex, FFTW_BACKWARD])



!MANY_VARIANT
#else



#ifndef MPI
        if (D%nthreads>1) call PoisFFT_InitThreads(1)
#endif
        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)

        call data_allocate_real(D%Solvers1D(1))

        D%Solvers1D(1)%forward = PoisFFT_Plan1D_Create(D%Solvers1D(1), [FFT_RealEven10])
        D%Solvers1D(1)%backward = PoisFFT_Plan1D_Create(D%Solvers1D(1), [FFT_RealEven01])

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          D%Solvers1D(i)%forward%planowner = .false.
          D%Solvers1D(i)%backward%planowner = .false.
          call data_allocate_real(D%Solvers1D(i))
        end do

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        call data_allocate_complex(D%Solvers2D(1))

        D%Solvers2D(1)%forward = PoisFFT_Plan2D_Create(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%backward = PoisFFT_Plan2D_Create(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])

        do i=2,D%nthreads
#ifdef MPI
          D%Solvers2D(i) = PoisFFT_Solver2D_From3D(D,3)
#else
          D%Solvers2D(i) = D%Solvers2D(1)
          D%Solvers2D(i)%forward%planowner = .false.
          D%Solvers2D(i)%backward%planowner = .false.
#endif
          call data_allocate_complex(D%Solvers2D(i))
#ifdef MPI
          D%Solvers2D(i)%forward = PoisFFT_Plan2D_Create(D%Solvers2D(i), [FFT_Complex, FFTW_FORWARD])
          D%Solvers2D(i)%backward = PoisFFT_Plan2D_Create(D%Solvers2D(i), [FFT_Complex, FFTW_BACKWARD])
#endif
        end do
        
        
        
#ifdef MPI
        call MPI_Comm_rank(D%Solvers1D(1)%mpi%comm, D%mpi%rank, ie)
        call MPI_Comm_size(D%Solvers1D(1)%mpi%comm, D%mpi%np, ie)
        allocate(D%mpi%snxs(D%mpi%np))
        allocate(D%mpi%snzs(D%mpi%np))
        allocate(D%mpi%rnxs(D%mpi%np))
        allocate(D%mpi%rnzs(D%mpi%np))
        allocate(D%mpi%sdispls(D%mpi%np))
        allocate(D%mpi%scounts(D%mpi%np))
        allocate(D%mpi%rdispls(D%mpi%np))
        allocate(D%mpi%rcounts(D%mpi%np))
        
        D%mpi%snxs(1:D%mpi%np-1) = (D%nx / D%mpi%np)
        D%mpi%snxs(D%mpi%np) = D%nx - sum(D%mpi%snxs(1:D%mpi%np-1))
        D%mpi%snzs = D%nz
        D%mpi%scounts = D%mpi%snxs*D%ny*D%nz
        D%mpi%sdispls(1) = 0
        do i = 2, D%mpi%np
          D%mpi%sdispls(i) = D%mpi%sdispls(i-1) + D%mpi%scounts(i-1)
        end do

        call MPI_AllToAll(D%mpi%snxs, 1, MPI_INTEGER, &
                          D%mpi%rnxs, 1, MPI_INTEGER, &
                          D%Solvers1D(1)%mpi%comm, ie)
        if (.not.all(D%mpi%rnxs(2:D%mpi%np)==D%mpi%rnxs(1))) stop "????"
        call MPI_AllToAll(D%mpi%snzs, 1, MPI_INTEGER, &
                          D%mpi%rnzs, 1, MPI_INTEGER, &
                          D%Solvers1D(1)%mpi%comm, ie)
        call MPI_AllToAll(D%mpi%scounts, 1, MPI_INTEGER, &
                          D%mpi%rcounts, 1, MPI_INTEGER, &
                          D%Solvers1D(1)%mpi%comm, ie)

        call MPI_AllToAll(D%mpi%sdispls, 1, MPI_INTEGER, &
                          D%mpi%rdispls, 1, MPI_INTEGER, &
                          D%Solvers1D(1)%mpi%comm, ie)

        D%mpi%rdispls(1) = 0
        do i = 2, D%mpi%np
          D%mpi%rdispls(i) = D%mpi%rdispls(i-1) + D%mpi%rcounts(i-1)
        end do
        !step1 local transpose
        allocate(D%mpi%tmp1(1:D%nz,1:D%ny,1:D%nx))
        allocate(D%mpi%tmp2(0:sum(D%mpi%rcounts)-1))
        allocate(D%mpi%rwork(sum(D%mpi%rnzs), D%ny, D%mpi%rnxs(1)))
#endif



!MANY_VARIANT
#endif

        forall(i=1:D%nx)  D%denomx(i) = 2 * (cos((i-1+D%offx)*dkx)-1.0_RP) / D%dx**2
        forall(j=1:D%ny)  D%denomy(j) = 2 * (cos((j-1+D%offy)*dky)-1.0_RP) / D%dy**2
        forall(k=1:D%nz)  D%denomz(k) = 2 * (cos((k-1+D%offz)*pi/D%gnz)-1.0_RP) / D%dz**2
      endif
    end subroutine PoisFFT_Solver3D_Init























    function PoisFFT_Solver2D_New(nxyz,dxyz,BCs,gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver2D) :: D

      integer, intent(in)   :: nxyz(2)
      real(RP), intent(in)  :: dxyz(2)
      integer, intent(in)   :: bcs(4)
      integer, intent(in), optional :: gnxyz(2)
      integer, intent(in), optional :: offs(2)
      integer, intent(in), optional  :: mpi_comm
      integer, intent(in), optional :: nthreads

      D%nxyz = nxyz

      D%dx = dxyz(1)
      D%dy = dxyz(2)

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

#ifdef MPI      
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
      else
        stop "No PFFT comm present in PoisFFT_Solver2D_New."
      end if
#endif

      !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver2D_New



    subroutine PoisFFT_Solver2D_Execute(D,Phi,RHS)
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
!
!       else if (all(D%BCs==PoisFFT_DirichletStag)) then
!
!         call PoisFFT_Solver2D_FullDirichletStag(D,&
!                  Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
!                      ngPhi(2)+1:ngPhi(2)+D%ny),&
!                  RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
!                      ngRHS(2)+1:ngRHS(2)+D%ny))

      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call PoisFFT_Solver2D_FullNeumannStag(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny))
      endif

    end subroutine PoisFFT_Solver2D_Execute





    subroutine PoisFFT_Solver2D_FullPeriodic(D, Phi, RHS)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:)
      real(RP), intent(in)  :: RHS(:,:)
      integer i,j

      ! Forward FFT of RHS
      D%cwork = cmplx(RHS,0._RP,CP)


#ifdef MPI
      call Execute_MPI(D%forward)
#else
      call Execute(D%forward, D%cwork)
#endif


      if (D%offx==0.and.D%offy==0) D%cwork(1,1) = 0

      do j = 1,D%ny
        do i = max(3-j-D%offx-D%offy,1),D%nx
          D%cwork(i,j) = D%cwork(i,j) / (D%denomx(i) + D%denomy(j))
        end do
      end do

#ifdef MPI
      call Execute_MPI(D%backward)
#else
      call Execute(D%backward, D%cwork)
#endif

      Phi = real(D%cwork,RP) / (D%gcnt)

    end subroutine PoisFFT_Solver2D_FullPeriodic



    subroutine PoisFFT_Solver2D_FullNeumannStag(D, Phi, RHS)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:)
      real(RP), intent(in)  :: RHS(:,:)
      integer i,j

      ! Forward FFT of RHS
      D%rwork = RHS


      call Execute(D%forward, D%rwork)


      D%rwork(1,1) = 0

      do j = 1,D%ny
        do i = max(3-j,1),D%nx
          D%rwork(i,j) = D%rwork(i,j) / (D%denomx(i) + D%denomy(j))
        end do
      end do

      call Execute(D%backward, D%rwork)

      Phi = D%rwork/(4*D%nx*D%ny)

    end subroutine PoisFFT_Solver2D_FullNeumannStag



    subroutine Poisfft_Solver2D_Init(D)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP) :: dkx, dky, dkx_h, dky_h
      integer :: i, j

      dkx_h = pi/(D%gnx)
      dky_h = pi/(D%gny)

      dkx = 2 * dkx_h
      dky = 2 * dky_h
        
      allocate(D%denomx(D%nx))
      allocate(D%denomy(D%ny))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call data_allocate_complex(D)

        D%forward = PoisFFT_Plan2D_Create(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan2D_Create(D, [FFT_Complex, FFTW_BACKWARD])

        forall(i=1:D%nx)  D%denomx(i) = 2 * (cos((i-1+D%offx)*dkx)-1.0_RP) / D%dx**2
        forall(j=1:D%ny)  D%denomy(j) = 2 * (cos((j-1+D%offy)*dky)-1.0_RP) / D%dy**2

      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan2D_Create(D, [(FFT_RealOdd10, i=1,2)])
        D%backward = PoisFFT_Plan2D_Create(D, [(FFT_RealOdd01, i=1,2)])

      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan2D_Create(D, [(FFT_RealEven10, i=1,2)])
        D%backward = PoisFFT_Plan2D_Create(D, [(FFT_RealEven01, i=1,2)])

        forall(i=1:D%nx)  D%denomx(i) = 2 * (cos((i-1)*dkx_h)-1.0_RP) / D%dx**2
        forall(j=1:D%ny)  D%denomy(j) = 2 * (cos((j-1)*dky_h)-1.0_RP) / D%dy**2

      endif
    end subroutine PoisFFT_Solver2D_Init






















    function PoisFFT_Solver1D_New(nxyz,dxyz,BCs,gnxyz,offs,mpi_comm,nthreads) result(D)
      type(PoisFFT_Solver1D) :: D

      integer, intent(in)   :: nxyz(1)
      real(RP), intent(in)  :: dxyz(1)
      integer, intent(in)   :: bcs(2)
      integer, intent(in), optional :: gnxyz(1)
      integer, intent(in), optional :: offs(1)
      integer, intent(in), optional :: mpi_comm
      integer, intent(in), optional :: nthreads
#ifdef MPI
      integer :: ie, dims
#endif
      D%nxyz = nxyz
      
      D%dx = dxyz(1)

      D%nx = nxyz(1)

      if (present(gnxyz)) then
        D%gnx = gnxyz(1)
      else
        D%gnx = D%nx
      end if

      if (present(offs)) then
        D%offx = offs(1)
      end if

      D%gcnt = int(D%gnx, kind(D%gcnt))

      D%cnt = D%nx

#ifdef MPI
      if (present(mpi_comm)) then
        D%mpi%comm = mpi_comm
      else
        stop "No PFFT comm present in PoisFFT_Solver1D_New."
      end if
      
      
      call MPI_Cartdim_get(D%mpi%comm, dims, ie)
      if (ie/=0) stop "Error executing MPI_Cartdim_get."
      D%mpi_transpose_needed = dims > 0
#endif

      D%BCs = BCs

      !create fftw plans and allocate working array
      call Init(D)
    end function PoisFFT_Solver1D_New



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
      integer i
#ifdef MPI
      interface
        subroutine MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
                              RECVTYPE, ROOT, COMM, IERROR)
           INTEGER    SENDBUF, RECVBUF(*)
           INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT
           INTEGER    COMM, IERROR
        end subroutine
        subroutine MPI_GATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS, &
                               DISPLS, RECVTYPE, ROOT, COMM, IERROR)
          import
          real(RP)    SENDBUF(*), RECVBUF(*)
          INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*)
          INTEGER    RECVTYPE, ROOT, COMM, IERROR
        end subroutine
        subroutine MPI_SCATTERV(SENDBUF, SENDCOUNTS, DISPLS, SENDTYPE, RECVBUF, &
                                RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
          import
          real(RP)    SENDBUF(*), RECVBUF(*)
          INTEGER    SENDCOUNTS(*), DISPLS(*), SENDTYPE
          INTEGER    RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
        end subroutine
      end interface
      integer, allocatable :: displs(:), counts(:)
      real(RP),allocatable :: tmp(:)
      integer :: ie

      if (D%mpi_transpose_needed) then
        call MPI_Comm_rank(D%mpi%comm, D%mpi%rank, ie)
        call MPI_Comm_size(D%mpi%comm, D%mpi%np, ie)
        allocate(displs(D%mpi%np))
        allocate(counts(D%mpi%np))
        if (D%mpi%rank==0) then
          allocate(tmp(1:D%gnx))
        else
          allocate(tmp(1))
        end if

        call MPI_Gather(D%offx, 1, MPI_INTEGER, &
                        displs, 1, MPI_INTEGER, 0, &
                        D%mpi%comm, ie)

        call MPI_Gather(D%nx, 1, MPI_INTEGER, &
                         counts, 1, MPI_INTEGER, 0, &
                         D%mpi%comm, ie)

        call MPI_Gatherv(RHS, D%nx, MPI_RP, &
                         tmp, counts, displs, MPI_RP, 0, &
                         D%mpi%comm, ie)
        if (ie/=0) stop "Error in MPI_Gatherv!"
        
        if (D%mpi%rank==0) D%cwork = cmplx(tmp,0._RP,CP)
      
        if (D%mpi%rank==0) then
          call Execute(D%forward, D%cwork)
        
          if (D%offx==0) D%cwork(1) = 0
          
          forall(i=2:D%gnx) &
            D%cwork(i) = D%cwork(i) / D%denomx(i)

          call Execute(D%backward, D%cwork)
        end if
      
        tmp = real(D%cwork,RP) / (D%nx)
        
        call MPI_Scatterv(tmp, counts, displs, MPI_RP, &
                          Phi, D%nx, MPI_RP, 0, &
                          D%mpi%comm, ie)
      end if
#else
      ! Forward FFT of RHS
      D%cwork = cmplx(RHS,0._RP,CP)

      call Execute(D%forward, D%cwork)

      forall(i=2:D%nx) &
        D%cwork(i) = D%cwork(i) / D%denomx(i)

      call Execute(D%backward, D%cwork)
      
      Phi = real(D%cwork,RP) / (D%nx)
#endif

    end subroutine PoisFFT_Solver1D_FullPeriodic





    subroutine PoisFFT_Solver1D_FullDirichletStag(D, Phi, RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      integer i

      ! Forward FFT of RHS
      D%rwork = RHS

      call Execute(D%forward, D%rwork)

      forall(i=1:D%nx) &
        D%rwork(i) = D%rwork(i) / D%denomx(i)

      call Execute(D%backward, D%rwork)

      Phi = D%rwork/(2*D%nx)

    end subroutine PoisFFT_Solver1D_FullDirichletStag




    subroutine PoisFFT_Solver1D_FullNeumannStag(D, Phi, RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      integer i

      ! Forward FFT of RHS
      D%rwork = RHS

      call Execute(D%forward, D%rwork)

      D%rwork = D%rwork/(D%nx)

      D%rwork(1) = 0
      forall(i=2:D%nx) &
        D%rwork(i) = D%rwork(i) / D%denomx(i)

      call Execute(D%backward, D%rwork)

      Phi = D%rwork/2

    end subroutine PoisFFT_Solver1D_FullNeumannStag



    subroutine Poisfft_Solver1D_Init(D)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP) :: dkx, dkx_h
      integer :: i

      dkx_h = pi/(D%gnx)

      dkx = 2 * dkx_h
      
      allocate(D%denomx(D%gnx))

      if (all(D%BCs==PoisFFT_Periodic)) then

        call data_allocate_complex(D)

        D%forward = PoisFFT_Plan1D_Create(D, [FFT_Complex, FFTW_FORWARD])
        D%backward = PoisFFT_Plan1D_Create(D, [FFT_Complex, FFTW_BACKWARD])
       
        forall(i=1:D%gnx) D%denomx(i) = 2 * (cos((i-1)*dkx)-1.0_RP) / D%dx**2

      else if (all(D%BCs==PoisFFT_DirichletStag)) then

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan1D_Create(D, [FFT_RealOdd10])
        D%backward = PoisFFT_Plan1D_Create(D, [FFT_RealOdd01])

        forall(i=1:D%gnx) D%denomx(i) = 2 * (cos((i)*dkx_h)-1.0_RP) / D%dx**2
        
      else if (all(D%BCs==PoisFFT_NeumannStag)) then

        call data_allocate_real(D)

        D%forward = PoisFFT_Plan1D_Create(D, [FFT_RealEven10])
        D%backward = PoisFFT_Plan1D_Create(D, [FFT_RealEven01])

        dkx_h = pi/(D%nx)
        
        forall(i=1:D%gnx) D%denomx(i) = 2 * (cos((i-1)*dkx_h)-1.0_RP) / D%dx**2
        
      endif
    end subroutine PoisFFT_Solver1D_Init

#undef RP
#undef CP
#undef _RP
#undef _CP
#undef MPI_RP
