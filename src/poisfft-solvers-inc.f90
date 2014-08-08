


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
      
        tmp = real(D%cwork,RP) / D%norm_factor
        
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
      
      Phi = real(D%cwork,RP) / D%norm_factor
#endif

    end subroutine PoisFFT_Solver1D_FullPeriodic





    subroutine PoisFFT_Solver1D_FullDirichlet(D, Phi, RHS)
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

      Phi = D%rwork / D%norm_factor

    end subroutine PoisFFT_Solver1D_FullDirichlet





    subroutine PoisFFT_Solver1D_FullNeumann(D, Phi, RHS)
      type(PoisFFT_Solver1D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:)
      real(RP), intent(in)  :: RHS(:)
      integer i

      ! Forward FFT of RHS
      D%rwork = RHS

      call Execute(D%forward, D%rwork)

      D%rwork(1) = 0
      forall(i=2:D%nx) &
        D%rwork(i) = D%rwork(i) / D%denomx(i)

      call Execute(D%backward, D%rwork)

      Phi = D%rwork / D%norm_factor

    end subroutine PoisFFT_Solver1D_FullNeumann









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

      Phi = real(D%cwork,RP) / D%norm_factor

    end subroutine PoisFFT_Solver2D_FullPeriodic



    subroutine PoisFFT_Solver2D_FullDirichlet(D, Phi, RHS)
      type(PoisFFT_Solver2D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:)
      real(RP), intent(in)  :: RHS(:,:)
      integer i,j

      ! Forward FFT of RHS
      D%rwork = RHS


      call Execute(D%forward, D%rwork)


      do j = 1,D%ny
        do i = 1,D%nx
          D%rwork(i,j) = D%rwork(i,j) / (D%denomx(i) + D%denomy(j))
        end do
      end do

      call Execute(D%backward, D%rwork)

      Phi = D%rwork / D%norm_factor

    end subroutine PoisFFT_Solver2D_FullDirichlet




    subroutine PoisFFT_Solver2D_FullNeumann(D, Phi, RHS)
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

      Phi = D%rwork / D%norm_factor
      
    end subroutine PoisFFT_Solver2D_FullNeumann




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
       Phi = real(D%cwork,RP) / D%norm_factor
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullPeriodic




    subroutine PoisFFT_Solver3D_FullDirichlet(D, Phi, RHS)
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
      Phi = D%rwork / D%norm_factor
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullDirichlet





    subroutine PoisFFT_Solver3D_FullNeumann(D, Phi, RHS)
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
      Phi = D%rwork / D%norm_factor
      !$omp end workshare

      !$omp end parallel

    end subroutine PoisFFT_Solver3D_FullNeumann




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
!               D%Solver1D%rwork(k,j,i) = real(D%Solver2D%cwork(i,j,k),RP) / D%norm_factor
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

        Phi(:,:,k) = real(D%Solvers2D(1)%cwork,RP) / D%norm_factor


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
!             D%Solver1D%rwork(k,j,i) = real(D%Solver2D%cwork(i,j,k),RP) / D%norm_factor
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

        Phi(:,:,k) = real(D%Solvers2D(tid)%cwork,RP) / D%norm_factor


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




