


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
!               mpi3%tmp1(k,j,i) = RHS(i,j,k)
!             end do
!           end do
!         end do
! 
!         
!         !step2 exchange
!         call MPI_AllToAllV(mpi3%tmp1, mpi3%scounts, mpi3%sdispls, MPI_RP, &
!                            mpi3%tmp2, mpi3%rcounts, mpi3%rdispls, MPI_RP, &
!                            D%Solver1mpi3%comm, ie)
!         !step3 local reordering of blocks
! 
!         do l = 1, mpi3%np
!           do k = 0, mpi3%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, mpi3%rnzs(l)-1
!                 D%Solver1D%rwork(i+sum(mpi3%rnzs(1:l-1))+1, j+1, k+1) = &
!                   mpi3%tmp2(i + j*mpi3%rnzs(l) + k*(D%ny*mpi3%rnzs(l)) + mpi3%rdispls(l))
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
!         do l = 1, mpi3%np
!           do k = 0, mpi3%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, mpi3%rnzs(l)-1
!                 mpi3%tmp2(i + j*mpi3%rnzs(l) + k*(D%ny*mpi3%rnzs(l)) + mpi3%rdispls(l)) = &
!                   D%Solver1D%rwork(i+sum(mpi3%rnzs(1:l-1))+1, j+1, k+1)
!               end do
!             end do
!           end do
!         end do
! 
!       
!         !step2' exchange
!         call MPI_AllToAllV(mpi3%tmp2, mpi3%rcounts, mpi3%rdispls, MPI_RP, &
!                            mpi3%tmp1, mpi3%scounts, mpi3%sdispls, MPI_RP, &
!                            D%Solver1mpi3%comm, ie)
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               D%Solver2D%cwork(i,j,k) = mpi3%tmp1(k,j,i)
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
!               mpi3%tmp1(k,j,i) = real(D%Solver2D%cwork(i,j,k),RP) / (2 * D%gcnt)
!             end do
!           end do
!         end do
! 
!         !step2 exchange
!         call MPI_AllToAllV(mpi3%tmp1, mpi3%scounts, mpi3%sdispls, MPI_RP, &
!                            mpi3%tmp2, mpi3%rcounts, mpi3%rdispls, MPI_RP, &
!                            D%Solver1mpi3%comm, ie)
! 
!         !step3 local reordering of blocks
!         do l = 1, mpi3%np
!           do k = 0, mpi3%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, mpi3%rnzs(l)-1
!                 D%Solver1D%rwork(i+sum(mpi3%rnzs(1:l-1))+1, j+1, k+1) = &
!                   mpi3%tmp2(i + j*mpi3%rnzs(l) + k*(D%ny*mpi3%rnzs(l)) + mpi3%rdispls(l))
!               end do
!             end do
!           end do
!         end do
!       
!         ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
!         call Execute(D%Solver1D%backward, D%Solver1D%rwork)
! 
!         !step3' local reordering of blocks
!         do l = 1, mpi3%np
!           do k = 0, mpi3%rnxs(1)-1
!             do j = 0, D%ny-1
!               do i = 0, mpi3%rnzs(l)-1
!                 mpi3%tmp2(i + j*mpi3%rnzs(l) + k*(D%ny*mpi3%rnzs(l)) + mpi3%rdispls(l)) = &
!                   D%Solver1D%rwork(i+sum(mpi3%rnzs(1:l-1))+1, j+1, k+1)
!               end do
!             end do
!           end do
!         end do
!       
!         !step2' exchange
!         call MPI_AllToAllV(mpi3%tmp2, mpi3%rcounts, mpi3%rdispls, MPI_RP, &
!                            mpi3%tmp1, mpi3%scounts, mpi3%sdispls, MPI_RP, &
!                            D%Solver1mpi3%comm, ie)
! 
!         do k = 1, D%nz
!           do j = 1, D%ny
!             do i = 1, D%nx
!               Phi(i,j,k) = mpi3%tmp1(k,j,i)
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

      
      if (D%Solvers1D(3)%mpi_transpose_needed) then
      
        call transform_1d_real_z(D%Solvers1D(3), Phi, RHS, forward=.true., use_rhs=.true.)
      
      else
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(3)%rwork = RHS(i,j,:)


            call Execute(D%Solvers1D(3)%forward,D%Solvers1D(3)%rwork)


            Phi(i,j,:) = D%Solvers1D(3)%rwork
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


      if (D%Solvers1D(3)%mpi_transpose_needed) then

        call transform_1d_real_z(D%Solvers1D(3), Phi, RHS, forward=.false., use_rhs=.false.)
        
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



#ifdef MPI
    subroutine PoisFFT_Solver3D_PNsNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout), target :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k

      if (D%Solvers1D(3)%mpi_transpose_needed) then
      
        call transform_1d_real_z(D%Solvers1D(3), Phi, RHS, forward=.true., use_rhs=.true.)
        
      else
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(3)%rwork = RHS(i,j,:)

            call Execute(D%Solvers1D(3)%forward,D%Solvers1D(3)%rwork)

            Phi(i,j,:) = D%Solvers1D(3)%rwork
          end do
        end do
      end if
      

      if (D%Solvers1D(2)%mpi_transpose_needed) then
      
        call transform_1d_real_y(D%Solvers1D(2), Phi, RHS, forward=.true., use_rhs=.false.)
        
      else
        do k = 1, D%nz
          do i = 1, D%nx
            D%Solvers1D(2)%rwork = Phi(i,:,k)

            call Execute(D%Solvers1D(2)%forward, D%Solvers1D(2)%rwork)
            
            Phi(i,:,k) = D%Solvers1D(2)%rwork
          end do
        end do
      end if

      do k = 1, D%nz
        do j = 1, D%ny
        
          D%Solvers1D(1)%cwork = cmplx(Phi(:,j,k),0._RP,CP)

          call Execute(D%Solvers1D(1)%forward, D%Solvers1D(1)%cwork)


          do i = max(4-j-k-D%offy-D%offz,1),D%nx
            D%Solvers1D(1)%cwork(i) = D%Solvers1D(1)%cwork(i) &
                                             / (D%denomx(i) + D%denomy(j) + D%denomz(k))
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          if (D%offy==0.and.D%offz==0.and.j==1.and.k==1) D%Solvers1D(1)%cwork(1) = 0

          call Execute(D%Solvers1D(1)%backward, D%Solvers1D(1)%cwork)

          Phi(:,j,k) = real(D%Solvers1D(1)%cwork,RP) / D%norm_factor

        end do
      end do

      if (D%Solvers1D(2)%mpi_transpose_needed) then

        call transform_1d_real_y(D%Solvers1D(2), Phi, RHS, forward=.false., use_rhs=.false.)

      else
        do k = 1, D%nz
          do i = 1, D%nx
            D%Solvers1D(2)%rwork = Phi(i,:,k)

            call Execute(D%Solvers1D(2)%backward, D%Solvers1D(2)%rwork)
            
            Phi(i,:,k) = D%Solvers1D(2)%rwork
          end do
        end do
      end if


      if (D%Solvers1D(3)%mpi_transpose_needed) then

        call transform_1d_real_z(D%Solvers1D(3), Phi, RHS, forward=.false., use_rhs=.false.)

      else
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(3)%rwork = Phi(i,j,:)


            call Execute(D%Solvers1D(3)%backward,D%Solvers1D(3)%rwork)


            Phi(i,j,:) = D%Solvers1D(3)%rwork
          end do
        end do
      end if
    end subroutine PoisFFT_Solver3D_PNsNs
    !MPI
#else
    !threads
    subroutine PoisFFT_Solver3D_PNsNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k
      integer tid      !thread id

      !$omp parallel private(tid,i,j,k)
      tid = 1
      !$ tid = omp_get_thread_num()+1

      ! Forward FFT of RHS in x dimension
      !$omp do
      do i = 1, D%nx
        D%Solvers2D(tid)%rwork = RHS(i,:,:)

        call Execute(D%Solvers2D(tid)%forward,D%Solvers2D(tid)%rwork)

        Phi(i,:,:) = D%Solvers2D(tid)%rwork
      end do
      !$omp end do

      !$omp do
      do k = 1, D%nz
        do j = 1, D%ny
        
          D%Solvers1D(tid)%cwork = cmplx(Phi(:,j,k),0._RP,CP)

          call Execute(D%Solvers1D(tid)%forward, D%Solvers1D(tid)%cwork)


          do i = max(4-j-k-D%offy-D%offz,1),D%nx
            D%Solvers1D(tid)%cwork(i) = D%Solvers1D(tid)%cwork(i) &
                                             / (D%denomx(i) + D%denomy(j) + D%denomz(k))
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          if (D%offy==0.and.D%offz==0.and.j==1.and.k==1) D%Solvers1D(tid)%cwork(1) = 0

          call Execute(D%Solvers1D(tid)%backward, D%Solvers1D(tid)%cwork)

          Phi(:,j,k) = real(D%Solvers1D(tid)%cwork,RP) / D%norm_factor

        end do
      end do
      !$omp end do

      !$omp do
      do i = 1, D%nx
        D%Solvers2D(tid)%rwork = Phi(i,:,:)

        call Execute(D%Solvers2D(tid)%backward,D%Solvers2D(tid)%rwork)

        Phi(i,:,:) = D%Solvers2D(tid)%rwork
      end do
      !$omp end do
      !$omp end parallel

    end subroutine PoisFFT_Solver3D_PNsNs
#endif
    
    
#ifdef MPI    
    subroutine transform_1d_real_y(D1D, Phi, RHS, forward, use_rhs)
      type(PoisFFT_Solver1D), intent(inout), target :: D1D
      real(RP), intent(inout) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      logical, intent(in) :: forward, use_rhs
      integer :: nx, ny, nz
      integer :: i,j,k,l, ie
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

      nx = size(Phi, 1)
      ny = size(Phi, 2)
      nz = size(Phi, 3)

      associate(mpi => D1D%mpi)
  
        !step1 local transpose
        if (use_rhs) then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                mpi%tmp1(j,k,i) = RHS(i,j,k)
              end do
            end do
          end do
        else
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                mpi%tmp1(j,k,i) = Phi(i,j,k)
              end do
            end do
          end do
        end if
        
        !step2 exchange
        call MPI_AllToAllV(mpi%tmp1, mpi%scounts, mpi%sdispls, MPI_RP, &
                           mpi%tmp2, mpi%rcounts, mpi%rdispls, MPI_RP, &
                           mpi%comm, ie)
                           
        !step3 local reordering of blocks
        do l = 1, mpi%np
          do k = 0, mpi%rnxs(1)-1
            do j = 0, nz-1
              do i = 0, mpi%rnzs(l)-1
                mpi%rwork(i+sum(mpi%rnzs(1:l-1))+1, j+1, k+1) = &
                  mpi%tmp2(i + j*mpi%rnzs(l) + k*(nz*mpi%rnzs(l)) + mpi%rdispls(l))
              end do
            end do
          end do
        end do
     
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        if (forward) then
          do k=1,size(mpi%rwork,3)
            do j=1,size(mpi%rwork,2)
              call Execute(D1D%forward,mpi%rwork(:,j,k))
            end do
          end do
        else
          do k=1,size(mpi%rwork,3)
            do j=1,size(mpi%rwork,2)
              call Execute(D1D%backward,mpi%rwork(:,j,k))
            end do
          end do
        end if      

        !step3' local reordering of blocks
        do l = 1, mpi%np
          do k = 0, mpi%rnxs(1)-1
            do j = 0, nz-1
              do i = 0, mpi%rnzs(l)-1
                mpi%tmp2(i + j*mpi%rnzs(l) + k*(nz*mpi%rnzs(l)) + mpi%rdispls(l)) = &
                  mpi%rwork(i+sum(mpi%rnzs(1:l-1))+1, j+1, k+1)
              end do
            end do
          end do
        end do

      
        !step2' exchange
        call MPI_AllToAllV(mpi%tmp2, mpi%rcounts, mpi%rdispls, MPI_RP, &
                           mpi%tmp1, mpi%scounts, mpi%sdispls, MPI_RP, &
                           mpi%comm, ie)
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              Phi(i,j,k) = mpi%tmp1(j,k,i)
            end do
          end do
        end do
        
      end associate
    end subroutine transform_1d_real_y
    
    
    subroutine transform_1d_real_z(D1D, Phi, RHS, forward, use_rhs)
      type(PoisFFT_Solver1D), intent(inout), target :: D1D
      real(RP), intent(inout) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      logical, intent(in) :: forward, use_rhs
      integer :: nx, ny, nz
      integer :: i,j,k,l, ie
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

      nx = size(Phi, 1)
      ny = size(Phi, 2)
      nz = size(Phi, 3)

      associate(mpi => D1D%mpi)
  
        !step1 local transpose
        if (use_rhs) then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                mpi%tmp1(k,j,i) = RHS(i,j,k)
              end do
            end do
          end do
        else
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                mpi%tmp1(k,j,i) = Phi(i,j,k)
              end do
            end do
          end do
        end if
        
        !step2 exchange
        call MPI_AllToAllV(mpi%tmp1, mpi%scounts, mpi%sdispls, MPI_RP, &
                           mpi%tmp2, mpi%rcounts, mpi%rdispls, MPI_RP, &
                           mpi%comm, ie)
                           
        !step3 local reordering of blocks
        do l = 1, mpi%np
          do k = 0, mpi%rnxs(1)-1
            do j = 0, ny-1
              do i = 0, mpi%rnzs(l)-1
                mpi%rwork(i+sum(mpi%rnzs(1:l-1))+1, j+1, k+1) = &
                  mpi%tmp2(i + j*mpi%rnzs(l) + k*(ny*mpi%rnzs(l)) + mpi%rdispls(l))
              end do
            end do
          end do
        end do
     
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        if (forward) then
          do k=1,size(mpi%rwork,3)
            do j=1,size(mpi%rwork,2)
              call Execute(D1D%forward,mpi%rwork(:,j,k))
            end do
          end do
        else
          do k=1,size(mpi%rwork,3)
            do j=1,size(mpi%rwork,2)
              call Execute(D1D%backward,mpi%rwork(:,j,k))
            end do
          end do
        end if      

        !step3' local reordering of blocks
        do l = 1, mpi%np
          do k = 0, mpi%rnxs(1)-1
            do j = 0, ny-1
              do i = 0, mpi%rnzs(l)-1
                mpi%tmp2(i + j*mpi%rnzs(l) + k*(ny*mpi%rnzs(l)) + mpi%rdispls(l)) = &
                  mpi%rwork(i+sum(mpi%rnzs(1:l-1))+1, j+1, k+1)
              end do
            end do
          end do
        end do

      
        !step2' exchange
        call MPI_AllToAllV(mpi%tmp2, mpi%rcounts, mpi%rdispls, MPI_RP, &
                           mpi%tmp1, mpi%scounts, mpi%sdispls, MPI_RP, &
                           mpi%comm, ie)
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              Phi(i,j,k) = mpi%tmp1(k,j,i)
            end do
          end do
        end do
        
      end associate
    end subroutine transform_1d_real_z

#endif    
    