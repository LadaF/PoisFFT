


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
    subroutine PoisFFT_Solver3D_PPNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k
      integer tid

      
      if (D%Solvers1D(3)%mpi_transpose_needed) then
      
        call transform_1d_real_z(D%Solvers1D(3:), Phi, RHS, forward=.true., use_rhs=.true.)
      
      else
        !$omp parallel private(tid,i,j,k)
        tid  = 0
        !$ tid = omp_get_thread_num()
        !$omp do collapse(2)
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(3+tid)%rwork = RHS(i,j,:)


            call Execute(D%Solvers1D(3+tid)%forward,D%Solvers1D(3+tid)%rwork)


            Phi(i,j,:) = D%Solvers1D(3+tid)%rwork
          end do
        end do
        !$omp end parallel
      end if
     
      if (D%ny<D%gny) then

        do k=1,D%nz

          !$omp parallel workshare
          D%Solvers2D(1)%cwork(1:D%nx,1:D%ny) = cmplx(Phi(1:D%nx,1:D%ny,k),0._RP,CP)
          !$omp end parallel workshare

          call Execute_MPI(D%Solvers2D(1)%forward, D%Solvers2D(1)%cwork)

          if (k==1.and.D%offz==0) then

            !$omp parallel do
            do j = 1,D%ny
              do i = max(3-j-D%offx-D%offy,1),D%nx
                D%Solvers2D(1)%cwork(i,j) = D%Solvers2D(1)%cwork(i,j) / (D%denomx(i) + D%denomy(j))
              end do
            end do
            !$omp end parallel do

            !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
            ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
            if (D%offx==0.and.D%offy==0) D%Solvers2D(1)%cwork(1,1) = 0

          else

            !$omp parallel do collapse(2)
            do j=1,D%ny
              do i=1,D%nx
                D%Solvers2D(1)%cwork(i,j) = D%Solvers2D(1)%cwork(i,j)&
                                                / (D%denomx(i) + D%denomy(j) + D%denomz(k))
              end do
            end do
            !$omp end parallel do

          endif

          call Execute_MPI(D%Solvers2D(1)%backward, D%Solvers2D(1)%cwork)

          !$omp parallel workshare
          Phi(:,:,k) = real(D%Solvers2D(1)%cwork,RP) / D%norm_factor
          !$omp end parallel workshare

        end do

      else

        !$omp parallel private(tid,i,j,k)
        tid  = 1
        !$ tid = omp_get_thread_num()+1

        !$omp do
        do k=1,D%nz
          D%Solvers2D(tid)%cwork(1:D%nx,1:D%ny) = cmplx(Phi(1:D%nx,1:D%ny,k),0._RP,CP)


          call Execute(D%Solvers2D(tid)%forward, D%Solvers2D(tid)%cwork)

          if (k==1.and.D%offz==0) then
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
        !$omp end parallel

      end if


      if (D%Solvers1D(3)%mpi_transpose_needed) then

        call transform_1d_real_z(D%Solvers1D(3:), Phi, RHS, forward=.false., use_rhs=.false.)
        
      else
        !$omp parallel private(tid,i,j,k)
        tid  = 0
        !$ tid = omp_get_thread_num()
        !$omp do collapse(2)
        do j=1,D%ny
          do i=1,D%nx
            D%Solvers1D(3+tid)%rwork = Phi(i,j,:)


            call Execute(D%Solvers1D(3+tid)%backward,D%Solvers1D(3+tid)%rwork)


            Phi(i,j,:) = D%Solvers1D(3+tid)%rwork
          end do
        end do
        !$omp end parallel
      end if
    end subroutine PoisFFT_Solver3D_PPNs
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




    subroutine PoisFFT_Solver3D_PNsP(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k
      integer tid      !thread id

      !$omp parallel private(tid,i,j,k)
      tid = 1
      !$ tid = omp_get_thread_num()+1

      ! Forward FFT of RHS in y dimension
      !$omp do
      do k=1,D%nz
        do i=1,D%nx
          D%Solvers1D(tid)%rwork = RHS(i,:,k)


          call Execute(D%Solvers1D(tid)%forward,D%Solvers1D(tid)%rwork)


          Phi(i,:,k) = D%Solvers1D(tid)%rwork
        end do
      end do
      !$omp end do

      
      !$omp do
      do j=1,D%ny
        D%Solvers2D(tid)%cwork(1:D%nx,1:D%nz) = cmplx(Phi(1:D%nx,1:D%nz,k),0._RP,CP)


        call Execute(D%Solvers2D(tid)%forward, D%Solvers2D(tid)%cwork)

        if (j==1) then
          do k=2,D%nz
            do i=2,D%nx
              D%Solvers2D(tid)%cwork(i,k) = D%Solvers2D(tid)%cwork(i,k)&
                                              / (D%denomx(i) + D%denomz(k))
            end do
          end do
          do i=2,D%nx
              D%Solvers2D(tid)%cwork(i,1) = D%Solvers2D(tid)%cwork(i,1)&
                                              / (D%denomx(i))
          end do
          do k=2,D%nz
              D%Solvers2D(tid)%cwork(1,k) = D%Solvers2D(tid)%cwork(1,k)&
                                              / (D%denomz(k))
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          D%Solvers2D(tid)%cwork(1,1) = 0

        else

          do k=1,D%nz
            do i=1,D%nx
              D%Solvers2D(tid)%cwork(i,k) = D%Solvers2D(tid)%cwork(i,k)&
                                              / (D%denomx(i) + D%denomy(j) + D%denomz(k))
            end do
          end do

        endif

        call Execute(D%Solvers2D(tid)%backward, D%Solvers2D(tid)%cwork)

        Phi(:,:,k) = real(D%Solvers2D(tid)%cwork,RP) / D%norm_factor


      end do
      !$omp end do

      
      !$omp do
      do k=1,D%nz
        do i=1,D%nx
          D%Solvers1D(tid)%rwork = Phi(i,:,k)


          call Execute(D%Solvers1D(tid)%backward,D%Solvers1D(tid)%rwork)


          Phi(i,:,k) = D%Solvers1D(tid)%rwork
        end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine PoisFFT_Solver3D_PNsP

    
    

#ifdef MPI
    subroutine PoisFFT_Solver3D_PNsNs(D, Phi, RHS)
      type(PoisFFT_Solver3D), intent(inout), target :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer i,j,k
      integer tid

      if (D%Solvers1D(3)%mpi_transpose_needed) then
      
        call transform_1d_real_z(D%Solvers1D(3::3), Phi, RHS, forward=.true., use_rhs=.true.)
        
      else
        !$omp parallel private(tid,i,j,k)
        tid  = 0
        !$ tid = omp_get_thread_num()
        !$omp do
        do j = 1, D%ny
          do i = 1, D%nx
            D%Solvers1D(3+3*tid)%rwork = RHS(i,j,:)

            call Execute(D%Solvers1D(3+3*tid)%forward,D%Solvers1D(3+3*tid)%rwork)

            Phi(i,j,:) = D%Solvers1D(3+3*tid)%rwork
          end do
        end do
        !$omp end parallel
     end if
      

      if (D%Solvers1D(2)%mpi_transpose_needed) then
      
        call transform_1d_real_y(D%Solvers1D(2::3), Phi, RHS, forward=.true., use_rhs=.false.)
        
      else
        !$omp parallel private(tid,i,j,k)
        tid  = 0
        !$ tid = omp_get_thread_num()
        !$omp do
        do k = 1, D%nz
          do i = 1, D%nx
            D%Solvers1D(2+3*tid)%rwork = Phi(i,:,k)

            call Execute(D%Solvers1D(2+3*tid)%forward, D%Solvers1D(2+3*tid)%rwork)
            
            Phi(i,:,k) = D%Solvers1D(2+3*tid)%rwork
          end do
        end do
        !$omp end parallel
      end if

      !$omp parallel private(tid,i,j,k)
      tid  = 0
      !$ tid = omp_get_thread_num()
      !$omp do
      do k = 1, D%nz
        do j = 1, D%ny
        
          D%Solvers1D(1+3*tid)%cwork = cmplx(Phi(:,j,k),0._RP,CP)

          call Execute(D%Solvers1D(1+3*tid)%forward, D%Solvers1D(1+3*tid)%cwork)


          do i = max(4-j-k-D%offy-D%offz,1),D%nx
            D%Solvers1D(1+3*tid)%cwork(i) = D%Solvers1D(1+3*tid)%cwork(i) &
                                             / (D%denomx(i) + D%denomy(j) + D%denomz(k))
          end do
          !NOTE: if IEEE FPE exceptions are disabled all this is not necessary and
          ! the loop can be over all indexes because the infinity or NaN is changed to 0 below
          if (D%offy==0.and.D%offz==0.and.j==1.and.k==1) D%Solvers1D(1+3*tid)%cwork(1) = 0

          call Execute(D%Solvers1D(1+3*tid)%backward, D%Solvers1D(1+3*tid)%cwork)

          Phi(:,j,k) = real(D%Solvers1D(1+3*tid)%cwork,RP) / D%norm_factor

        end do
      end do
      !$omp end parallel

      if (D%Solvers1D(2)%mpi_transpose_needed) then

        call transform_1d_real_y(D%Solvers1D(2::3), Phi, RHS, forward=.false., use_rhs=.false.)

      else
        !$omp parallel private(tid,i,j,k)
        tid  = 0
        !$ tid = omp_get_thread_num()
        !$omp do
        do k = 1, D%nz
          do i = 1, D%nx
            D%Solvers1D(2+3*tid)%rwork = Phi(i,:,k)

            call Execute(D%Solvers1D(2+3*tid)%backward, D%Solvers1D(2+3*tid)%rwork)
            
            Phi(i,:,k) = D%Solvers1D(2+3*tid)%rwork
          end do
        end do
        !$omp end parallel
      end if


      if (D%Solvers1D(3)%mpi_transpose_needed) then

        call transform_1d_real_z(D%Solvers1D(3::3), Phi, RHS, forward=.false., use_rhs=.false.)

      else
       !$omp parallel private(tid,i,j,k)
        tid  = 0
        !$ tid = omp_get_thread_num()
        !$omp do
        do j = 1, D%ny
          do i = 1, D%nx
            D%Solvers1D(3+3*tid)%rwork = Phi(i,j,:)

            call Execute(D%Solvers1D(3+3*tid)%backward,D%Solvers1D(3+3*tid)%rwork)

            Phi(i,j,:) = D%Solvers1D(3+3*tid)%rwork
          end do
        end do
        !$omp end parallel
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
      type(PoisFFT_Solver1D), intent(inout), target :: D1D(:)
      real(RP), intent(inout) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      logical, intent(in) :: forward, use_rhs
      integer :: nx, ny, nz
      integer :: i,j,k,l, ie
      integer :: tid
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

#define m D1D(1)%mpi

        !$omp parallel private(i,j,k,tid)
        tid = 1
        !$ tid = omp_get_thread_num()+1
   
        !step1 local transpose
        if (use_rhs) then
          !$omp do
          do k = 1, nz
            do i = 1, nx
              do j = 1, ny
                m%tmp1(j,k,i) = RHS(i,j,k)
              end do
            end do
          end do
        else
          !$omp do
          do k = 1, nz
            do i = 1, nx
              do j = 1, ny
                m%tmp1(j,k,i) = Phi(i,j,k)
              end do
            end do
          end do
        end if
        
        !step2 exchange
        !$omp single
        call MPI_AllToAllV(m%tmp1, m%scounts, m%sdispls, MPI_RP, &
                           m%tmp2, m%rcounts, m%rdispls, MPI_RP, &
                           m%comm, ie)
        !$omp end single
                           
        !step3 local reordering of blocks
        !$omp do collapse(3)
        do l = 1, m%np
          do k = 0, m%rnxs(1)-1
            do j = 0, nz-1
              do i = 0, m%rnzs(l)-1
                m%rwork(i+m%sumrnzs(l)+1, j+1, k+1) = &
                  m%tmp2(i + j*m%rnzs(l) + k*(nz*m%rnzs(l)) + m%rdispls(l))
              end do
            end do
          end do
        end do
     
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        if (forward) then
          !$omp do
          do k=1,size(m%rwork,3)
            do j=1,size(m%rwork,2)
              !TODO: is this really alligned?
              call Execute(D1D(tid)%forward,m%rwork(:,j,k))
            end do
          end do
        else
          !$omp do
          do k=1,size(m%rwork,3)
            do j=1,size(m%rwork,2)
              call Execute(D1D(tid)%backward,m%rwork(:,j,k))
            end do
          end do
        end if      

        !step3' local reordering of blocks
        !$omp do collapse(3)
        do l = 1, m%np
          do k = 0, m%rnxs(1)-1
            do j = 0, nz-1
              do i = 0, m%rnzs(l)-1
                m%tmp2(i + j*m%rnzs(l) + k*(nz*m%rnzs(l)) + m%rdispls(l)) = &
                  m%rwork(i+m%sumrnzs(l)+1, j+1, k+1)
              end do
            end do
          end do
        end do

      
        !step2' exchange
        !$omp single
        call MPI_AllToAllV(m%tmp2, m%rcounts, m%rdispls, MPI_RP, &
                           m%tmp1, m%scounts, m%sdispls, MPI_RP, &
                           m%comm, ie)
        !$omp end single

        !$omp do
        do k = 1, nz
          do i = 1, nx
            do j = 1, ny
              Phi(i,j,k) = m%tmp1(j,k,i)
            end do
          end do
        end do

        !$omp end parallel
        
#undef mpi
    end subroutine transform_1d_real_y
    
    
    subroutine transform_1d_real_z(D1D, Phi, RHS, forward, use_rhs)
      type(PoisFFT_Solver1D), intent(inout), target :: D1D(:)
      real(RP), intent(inout) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      logical, intent(in) :: forward, use_rhs
      integer :: nx, ny, nz
      integer :: i,j,k,l, ie
      integer :: tid
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

#define m D1D(1)%mpi
  
        !$omp parallel private(i,j,k,tid)
        tid = 1
        !$ tid = omp_get_thread_num()+1

        !step1 local transpose
        if (use_rhs) then
         !$omp do
          do j = 1, ny
            do i = 1, nx
              do k = 1, nz
                m%tmp1(k,j,i) = RHS(i,j,k)
              end do
            end do
          end do
        else
         !$omp do
          do j = 1, ny
            do i = 1, nx
              do k = 1, nz
                m%tmp1(k,j,i) = Phi(i,j,k)
              end do
            end do
          end do
        end if
        
        !step2 exchange
        !$omp single
        call MPI_AllToAllV(m%tmp1, m%scounts, m%sdispls, MPI_RP, &
                           m%tmp2, m%rcounts, m%rdispls, MPI_RP, &
                           m%comm, ie)
        !$omp end single
                           
        !step3 local reordering of blocks
        !$omp do collapse(3)
        do l = 1, m%np
          do k = 0, m%rnxs(1)-1
            do j = 0, ny-1
              do i = 0, m%rnzs(l)-1
                m%rwork(i+m%sumrnzs(l)+1, j+1, k+1) = &
                  m%tmp2(i + j*m%rnzs(l) + k*(ny*m%rnzs(l)) + m%rdispls(l))
              end do
            end do
          end do
        end do
     
        ! Forward FFT of RHS in z dimension according to Wilhelmson, Ericksen, JCP 1977
        if (forward) then
          !$omp do
          do k=1,size(m%rwork,3)
            do j=1,size(m%rwork,2)
              call Execute(D1D(tid)%forward,m%rwork(:,j,k))
            end do
          end do
        else
          !$omp do
          do k=1,size(m%rwork,3)
            do j=1,size(m%rwork,2)
              call Execute(D1D(tid)%backward,m%rwork(:,j,k))
            end do
          end do
        end if      

        !step3' local reordering of blocks
        !$omp do collapse(3)
        do l = 1, m%np
          do k = 0, m%rnxs(1)-1
            do j = 0, ny-1
              do i = 0, m%rnzs(l)-1
                m%tmp2(i + j*m%rnzs(l) + k*(ny*m%rnzs(l)) + m%rdispls(l)) = &
                  m%rwork(i+m%sumrnzs(l)+1, j+1, k+1)
              end do
            end do
          end do
        end do

      
        !step2' exchange
        !$omp single
        call MPI_AllToAllV(m%tmp2, m%rcounts, m%rdispls, MPI_RP, &
                           m%tmp1, m%scounts, m%sdispls, MPI_RP, &
                           m%comm, ie)
        !$omp end single

        !$omp do
        do j = 1, ny
          do k = 1, nz
            do i = 1, nx
              Phi(i,j,k) = m%tmp1(k,j,i)
            end do
          end do
        end do
        
       !$omp end parallel
        
#undef mpi
    end subroutine transform_1d_real_z

#endif    
    