module PoisFFT
  use iso_c_binding
  !$ use omp_lib
  use FFT
  implicit none

  integer, parameter :: PoisFFT_Periodic = 0
  integer, parameter :: PoisFFT_Dirichlet = 1
  integer, parameter :: PoisFFT_Neumann = 2
  integer, parameter :: PoisFFT_DirichletStag = 3
  integer, parameter :: PoisFFT_NeumannStag = 4

  real(RP), parameter, private :: pi = 3.141592653589793238462_RP

  interface DeallocateData
    module procedure PoisFFT_Solver1D_DeallocateData
    module procedure PoisFFT_Solver2D_DeallocateData
    module procedure PoisFFT_Solver3D_DeallocateData
  end interface DeallocateData
  
  contains

    subroutine poisFFT_Solver3D_Execute(D,Phi,RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)

      integer   :: ngPhi(3), ngRHS(3)


      ngPhi = (ubound(Phi)-[D%nx,D%ny,D%nz])/2
      ngRHS = (ubound(RHS)-[D%nx,D%ny,D%nz])/2


      if (all(D%BCs==PoisFFT_Periodic)) then
write(*,*) "FullPeriodic",D%BCs      
        call PoisFFT_Solver3D_FullPeriodic(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))
                     
      else if (all(D%BCs(1:4)==PoisFFT_Periodic).and.all(D%BCs(5:6)==PoisFFT_NeumannStag)) then
write(*,*) "PPNs",D%BCs
        call PoisFFT_Solver3D_PPNs(D,&
                 Phi(ngPhi(1)+1:ngPhi(1)+D%nx,&
                     ngPhi(2)+1:ngPhi(2)+D%ny,&
                     ngPhi(3)+1:ngPhi(3)+D%nz),&
                 RHS(ngRHS(1)+1:ngRHS(1)+D%nx,&
                     ngRHS(2)+1:ngRHS(2)+D%ny,&
                     ngRHS(3)+1:ngRHS(3)+D%nz))

       endif

    end subroutine poisFFT_Solver3D_Execute


    
    function PoisFFT_Solver3D_New(nx,ny,nz,dx,dy,dz,BCs) result(D)
      type(PoisFFT_Solver3D) :: D

      integer, intent(in)   :: nx, ny, nz
      real(RP), intent(in)  :: dx, dy, dz
      integer, dimension(6) :: bcs

      D%dx = dx
      D%dy = dy
      D%dz = dz

      D%nx = nx
      D%ny = ny
      D%nz = nz

      D%cnt = D%nx * D%ny * D%nz

      D%BCs = BCs
write(*,*) "BCs",D%BCs
      !create fftw plans and allocate working array
      call poisfft_3d_init(D)
    end function PoisFFT_Solver3D_New


    function PoisFFT_Solver1D_From3D(D3D,direction) result(D)
      type(PoisFFT_Solver1D)             :: D
      type(PoisFFT_Solver3D), intent(in) :: D3D
      integer, intent(in)                :: direction

      if (direction==1) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%BCs = D3D%BCs(1:2)
      else if (direction==2) then
        D%dx = D3D%dy
        D%nx = D3D%ny
        D%BCs = D3D%BCs(3:4)
      else
        D%dx = D3D%dz
        D%nx = D3D%nz
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
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%BCs = D3D%BCs(3:6)
      else if (direction==2) then
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%dy = D3D%dz
        D%ny = D3D%nz
        D%BCs = D3D%BCs([1,2,5,6])
      else
        D%dx = D3D%dx
        D%nx = D3D%nx
        D%dy = D3D%dy
        D%ny = D3D%ny
        D%BCs = D3D%BCs(1:4)
      endif

      D%cnt = D%nx * D%ny
    end function PoisFFT_Solver2D_From3D


    subroutine PoisFFT_Solver1D_DeallocateData(D)
      type(PoisFFT_Solver1D), intent(inout) :: D

      call Destroy(D%fplan)
      call Destroy(D%bplan)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

    endsubroutine PoisFFT_Solver1D_DeallocateData



    subroutine PoisFFT_Solver2D_DeallocateData(D)
      type(PoisFFT_Solver2D), intent(inout) :: D

      call Destroy(D%fplan)
      call Destroy(D%bplan)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

    endsubroutine PoisFFT_Solver2D_DeallocateData



    subroutine PoisFFT_Solver3D_DeallocateData(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer :: i

      call Destroy(D%fplan)
      call Destroy(D%bplan)

      if (associated(D%rwork)) call data_deallocate(D%rwork)
      if (associated(D%cwork)) call data_deallocate(D%cwork)

      if (allocated(D%Solvers1D)) then
        do i=1,size(D%Solvers1D)
          call DeallocateData(D%Solvers1D(i))
        enddo
      endif

      if (allocated(D%Solvers2D)) then
        do i=1,size(D%Solvers2D)
          call DeallocateData(D%Solvers2D(i))
        enddo
      endif

      if (allocated(D%Solvers2D)) then
        do i=1,size(D%Solvers2D)
          call DeallocateData(D%Solvers2D(i))
        enddo
      endif
    endsubroutine PoisFFT_Solver3D_DeallocateData



    subroutine PoisFFT_Solver3D_FullPeriodic(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      real(RP) :: dx2, dy2, dz2
      real(RP) :: dkx, dky, dkz
      integer i,j,k

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2
      dz2 = 1._RP/D%dz**2

      dkx = 2._RP*pi/(D%nx)
      dky = 2._RP*pi/(D%ny)
      dkz = 2._RP*pi/(D%nz)

      ! Forward FFT of RHS
      forall(i=1:D%nx,j=1:D%ny,k=1:D%nz)&
        D%cwork(i,j,k) = cmplx(RHS(i,j,k),0._RP)


      call Execute(D%fplan, D%cwork)

      D%cwork(1,1,1) = 0
      forall(i=1:D%nx,j=1:D%ny,k=1:D%nz, i/=1.or.j/=1.or.k/=1) &
        D%cwork(i,j,k) = D%cwork(i,j,k) / (2.0_RP * ((cos((i-1)*dkx)-1.0_RP)*dx2+&
                                                   (cos((j-1)*dky)-1.0_RP)*dy2+&
                                                   (cos((k-1)*dkz)-1.0_RP)*dz2))

      call Execute(D%bplan, D%cwork)

      Phi = real(D%cwork,RP)/(D%nx*D%ny*D%nz)

    end subroutine PoisFFT_Solver3D_FullPeriodic




    subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      real(RP) :: dx2, dy2, dz2
      real(RP) :: dkx, dky, dkz
      integer i,j,k
      integer tid      !thread id

      dx2 = 1._RP/D%dx**2
      dy2 = 1._RP/D%dy**2
      dz2 = 1._RP/D%dz**2

      dkx = 2._RP*pi/(D%nx)
      dky = 2._RP*pi/(D%ny)
      dkz = 2._RP*pi/(D%nz)

      !$omp parallel private(tid)
      tid = 1
      !$ tid = omp_get_thread_num()
      

      ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
      !$omp do
      do j=1,D%ny
        do i=1,D%nx
          D%Solvers1D(tid)%rwork = RHS(i,j,:)
          
          call Execute(D%Solvers1D(tid)%fplan,D%Solvers1D(tid)%rwork)
          
          Phi(i,j,:) = D%Solvers1D(tid)%rwork
        enddo
      enddo
      !$omp end do

      !now solve nz 2D problems
      !$omp do
      do k=1,D%nz
        forall(i=1:D%nx,j=1:D%ny)&
          D%Solvers2D(tid)%cwork(i,j) = cmplx(Phi(i,j,k),0._RP)
          
        call Execute(D%Solvers2D(tid)%fplan, D%Solvers2D(tid)%cwork)

        D%Solvers2D(tid)%cwork(1,1) = 0
        forall(i=1:D%nx,j=1:D%ny, i/=1.or.j/=1) &
          D%Solvers2D(tid)%cwork(i,j) = D%Solvers2D(tid)%cwork(i,j)&
                                            / (2.0_RP * ((cos((i-1)*dkx)-1.0_RP)*dx2+&
                                                        (cos((j-1)*dky)-1.0_RP)*dy2+&
                                                        (cos((k-1)*pi/D%nz)-1.0_RP)*dz2))

        call Execute(D%Solvers2D(tid)%bplan, D%Solvers2D(tid)%cwork)

        Phi(:,:,k) = real(D%Solvers2D(tid)%cwork,RP)/(D%nx*D%ny)

        if (k>1) then
          Phi(:,:,k) = 2 * Phi(:,:,k) / D%nz
        else
          Phi(:,:,k) = Phi(:,:,k) / D%nz
        endif
      enddo
      !$omp end do

      !$omp do
      do j=1,D%ny
        do i=1,D%nx
          D%Solvers1D(tid)%rwork = Phi(i,j,:)
          
          call Execute(D%Solvers1D(tid)%bplan,D%Solvers1D(tid)%rwork)
          
          Phi(i,j,:) = D%Solvers1D(tid)%rwork / 4._RP
        enddo
      enddo
      !$omp end do

      !$omp end parallel
    end subroutine PoisFFT_Solver3D_PPNs




    subroutine poisfft_3d_init(D)
      type(PoisFFT_Solver3D), intent(inout) :: D
      integer i
      
      D%nthreads = 1
      !$omp parallel
      !$omp critical
      !$ D%nthreads = omp_get_num_threads()
      !$omp end critical
      !$omp end parallel
write (*,*) "init",D%BCs      
      if (all(D%BCs==PoisFFT_Periodic)) then
      
        D%fplan = PoisFFT_Plan3D_QuickCreate(D, [FFT_Complex, FFTW_FORWARD])
        D%bplan = PoisFFT_Plan3D_QuickCreate(D, [FFT_Complex, FFTW_BACKWARD])
        
        call data_allocate_complex(D)
        
      else if (all(D%BCs(1:4)==PoisFFT_Periodic).and.all(D%BCs(5:6)==PoisFFT_NeumannStag)) then

        allocate(D%Solvers1D(D%nthreads))

        D%Solvers1D(1) = PoisFFT_Solver1D_From3D(D,3)
        
        D%Solvers1D(1)%fplan = PoisFFT_Plan1D_QuickCreate(D%Solvers1D(1), [FFT_RealEven10])
        D%Solvers1D(1)%bplan = PoisFFT_Plan1D_QuickCreate(D%Solvers1D(1), [FFT_RealEven01])
        call data_allocate_real(D%Solvers1D(1))

        do i=2,D%nthreads
          D%Solvers1D(i) = D%Solvers1D(1)
          call data_allocate_real(D%Solvers1D(i))
        enddo

        allocate(D%Solvers2D(D%nthreads))

        D%Solvers2D(1) = PoisFFT_Solver2D_From3D(D,3)

        D%Solvers2D(1)%fplan = PoisFFT_Plan2D_QuickCreate(D%Solvers2D(1), [FFT_Complex, FFTW_FORWARD])
        D%Solvers2D(1)%bplan = PoisFFT_Plan2D_QuickCreate(D%Solvers2D(1), [FFT_Complex, FFTW_BACKWARD])
        call data_allocate_complex(D%Solvers2D(1))

        do i=2,D%nthreads
          D%Solvers2D(i) = D%Solvers2D(1)
          call data_allocate_complex(D%Solvers2D(i))
        enddo

      endif
    end subroutine poisfft_3d_init

    
end module PoisFFT
