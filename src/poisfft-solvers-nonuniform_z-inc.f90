    subroutine PoisFFT_Solver3D_nonuniform_Z_PPNs(D, Phi, RHS)
      type(PoisFFT_Solver3D_nonuniform_Z), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer :: i, j, k
      integer :: tid      !thread id
      real(RP) :: lam
      
      tid = 1
! #ifdef MPI
!       do k = 1, D%nz
!           D%Solvers2D(tid)%cwork = RHS(:,:,k)
! 
!           call Execute_MPI(D%Solvers2D(tid)%forward)
! 
!           D%cwork(:,:,k) = D%Solvers2D(tid)%cwork
!       end do
! 
!       if (D%Solvers1D(3)%mpi_transpose_needed) then
!         call solve_tridiagonal_1d_complex_z(D%Solvers1D(3:), D, Phi)
!       else
!         !$omp do private(lam)
!         do j = 1, D%ny
!           do i = 1, D%nx
!             lam = D%denomx(i) + D%denomy(j) 
!             if (i==1.and.j==1) then
!               call solve_tridiag(1._RP, 0._RP, lam, D%cwork(i,j,:))
!             else          
!               call solve_tridiag(D%mat_b(1), D%mat_c(1), lam, D%cwork(i,j,:))
!             end if
!           end do
!         end do      
!         !$omp end do
!       endif
! 
!       do k = 1, D%nz
!           D%Solvers2D(tid)%cwork = D%cwork(:,:,k)
!           
!           call Execute_MPI(D%Solvers2D(tid)%backward)          
!           
!           Phi(:,:,k) = real(D%Solvers2D(tid)%cwork, kind=RP) / (D%norm_factor)
!       end do
! 
! 
! #else

      !$omp parallel private(tid,i,j,k)
      !$ tid = omp_get_thread_num()+1

      !$omp do
      do k = 1, D%nz
          D%Solvers2D(tid)%cwork = RHS(:,:,k)

          call Execute(D%Solvers2D(tid)%forward,D%Solvers2D(tid)%cwork)
          
          D%cwork(:,:,k) = D%Solvers2D(tid)%cwork
      end do
      !$omp end do

      !$omp do private(lam)
      do j = 1, D%ny
        do i = 1, D%nx
          lam = D%denomx(i) + D%denomy(j) 
          if (i==1.and.j==1) then
            call solve_tridiag(D%mat_b(1), 0._RP, lam, D%cwork(i,j,:))
          else          
            call solve_tridiag(D%mat_b(1), D%mat_c(1), lam, D%cwork(i,j,:))
          end if
        end do
      end do      
      !$omp end do

      !$omp do
      do k = 1, D%nz
          D%Solvers2D(tid)%cwork = D%cwork(:,:,k)
          
          call Execute(D%Solvers2D(tid)%backward,D%Solvers2D(tid)%cwork)          
          
          Phi(:,:,k) = real(D%Solvers2D(tid)%cwork, kind=RP) / (D%norm_factor)
      end do
      !$omp end do
      !$omp end parallel
! #endif     
    contains
    

      subroutine solve_tridiag(b1, c1, lambdas, x)
          real(RP), intent(in) :: b1, c1, lambdas
          complex(CP), intent(inout) :: x(:)
          real(RP) :: cp(D%gnz)
          real(RP) :: m
          integer :: i, n

          n = D%gnz

  ! initialize c-prime
          m = b1 + lambdas
          x(1) = x(1) / m
          if (n==1) return
          cp(1) = c1 / m
  ! solve for vectors c-prime and x-prime
          do i = 2, n-1
            m = D%mat_b(i) + lambdas - cp(i-1)*D%mat_a(i)
            cp(i) = D%mat_c(i)/m
            x(i) = (x(i)-x(i-1)*D%mat_a(i)) / m
          enddo
          x(n) = (x(n)-x(n-1)*D%mat_a(n)) / &
                 (D%mat_b(n) + lambdas - cp(n-1)*D%mat_a(n))
  ! solve for x from the vectors c-prime and x-prime
          do i = n-1, 1, -1
            x(i) = x(i) - cp(i)*x(i+1)
          end do

      end subroutine solve_tridiag 
    end subroutine PoisFFT_Solver3D_nonuniform_Z_PPNs

  
    subroutine PoisFFT_Solver3D_nonuniform_Z_FullNeumann(D, Phi, RHS)
      type(PoisFFT_Solver3D_nonuniform_Z), intent(inout) :: D
      real(RP), intent(out) :: Phi(:,:,:)
      real(RP), intent(in)  :: RHS(:,:,:)
      integer :: i, j, k
      integer :: tid      !thread id
      real(RP) :: lam
      
      tid = 1
#ifdef MPI
      do k = 1, D%nz
          D%Solvers2D(tid)%rwork = RHS(:,:,k)

          call Execute_MPI(D%Solvers2D(tid)%forward)
          
          Phi(:,:,k) = D%Solvers2D(tid)%rwork
      end do

      if (D%Solvers1D(3)%mpi_transpose_needed) then
        call solve_tridiagonal_1d_real_z(D%Solvers1D(3:), D, Phi)
      else
        !$omp parallel do private(lam) collapse(2)
        do j = 1, D%ny
          do i = 1, D%nx
            lam = D%denomx(i) + D%denomy(j)
            if (i+D%offx==1.and.j+D%offy==1) then
              call solve_tridiag(D%mat_b(1), 0._RP, lam, Phi(i,j,:))
            else          
              call solve_tridiag(D%mat_b(1), D%mat_c(1), lam, Phi(i,j,:))
            end if
          end do
        end do      
        !$omp end parallel do
      end if

      do k = 1, D%nz
          D%Solvers2D(tid)%rwork = Phi(:,:,k)
          
          call Execute_MPI(D%Solvers2D(tid)%backward)
          
          Phi(:,:,k) = real(D%Solvers2D(tid)%rwork, kind=RP) / (D%norm_factor)
      end do
      
      
#else

      !$omp parallel private(tid,i,j,k)
      !$ tid = omp_get_thread_num()+1

      !$omp do
      do k = 1, D%nz
          D%Solvers2D(tid)%rwork = RHS(:,:,k)

          call Execute(D%Solvers2D(tid)%forward,D%Solvers2D(tid)%rwork)
          
          Phi(:,:,k) = D%Solvers2D(tid)%rwork
      end do
      !$omp end do

      !$omp do private(lam) collapse(2)
      do j = 1, D%ny
        do i = 1, D%nx
          lam = D%denomx(i) + D%denomy(j)
          if (i==1.and.j==1) then
            call solve_tridiag(D%mat_b(1), 0._RP, lam, Phi(i,j,:))
          else          
            call solve_tridiag(D%mat_b(1), D%mat_c(1), lam, Phi(i,j,:))
          end if
        end do
      end do      
      !$omp end do

      !$omp do
      do k = 1, D%nz
          D%Solvers2D(tid)%rwork = Phi(:,:,k)
          
          call Execute(D%Solvers2D(tid)%backward,D%Solvers2D(tid)%rwork)
          
          Phi(:,:,k) = real(D%Solvers2D(tid)%rwork, kind=RP) / (D%norm_factor)
      end do
      !$omp end do
      !$omp end parallel
#endif
    contains
    

      subroutine solve_tridiag(b1, c1, lambdas, x)
          real(RP), intent(in) :: b1, c1, lambdas
          real(RP), intent(inout) :: x(:)
          real(RP) :: cp(D%gnz)
          real(RP) :: m
          integer :: i, n

          n = D%gnz

  ! initialize c-prime
          m = b1 + lambdas
          x(1) = x(1) / m
          if (n==1) return
          cp(1) = c1 / m
  ! solve for vectors c-prime and x-prime
          do i = 2, n-1
            m = D%mat_b(i) + lambdas - cp(i-1)*D%mat_a(i)
            cp(i) = D%mat_c(i)/m
            x(i) = (x(i)-x(i-1)*D%mat_a(i)) / m
          enddo
          x(n) = (x(n)-x(n-1)*D%mat_a(n)) / &
                 (D%mat_b(n) + lambdas - cp(n-1)*D%mat_a(n))
  ! solve for x from the vectors c-prime and x-prime
          do i = n-1, 1, -1
            x(i) = x(i) - cp(i)*x(i+1)
          end do

      end subroutine solve_tridiag 
    end subroutine PoisFFT_Solver3D_nonuniform_Z_FullNeumann


    
    


#ifdef MPI
    subroutine solve_tridiagonal_1d_real_z(D1D, D, Phi)
      type(PoisFFT_Solver1D), intent(inout), target :: D1D(:)
      type(PoisFFT_Solver3D_nonuniform_Z), intent(inout), target :: D
      real(RP), intent(inout) :: Phi(:,:,:)
      integer :: nx, ny, nz
      integer :: i,j,k,l, ie
      integer :: glob_i
      integer :: tid
      real(RP) :: lam
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
        !$omp do collapse(3)
        do j = 1, ny
          do i = 1, nx
            do k = 1, nz
              m%tmp1(k,j,i) = Phi(i,j,k)
            end do
          end do
        end do

        
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
     
        !$omp do collapse(2)
        do k=1,size(m%rwork,3)
          do j=1,size(m%rwork,2)
            glob_i = k + D1D(1)%mpi%rank * (D%nx / D1D(1)%mpi%np)
            lam = D%denomx(glob_i) + D%denomy(j)
            if (glob_i==1.and.j+D%offy==1) then
                call solve_tridiag(D%mat_b(1), 0._RP, lam, m%rwork(:,j,k))
            else          
                call solve_tridiag(D%mat_b(1), D%mat_c(1), lam, m%rwork(:,j,k))
            end if
          end do
        end do

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

        !$omp do collapse(3)
        do j = 1, ny
          do k = 1, nz
            do i = 1, nx
              Phi(i,j,k) = m%tmp1(k,j,i)
            end do
          end do
        end do
        
       !$omp end parallel
        
#undef m       
    contains
    
      subroutine solve_tridiag(b1, c1, lambdas, x)
          real(RP), intent(in) :: b1, c1, lambdas
          real(RP), intent(inout) :: x(:)
          real(RP) :: cp(D%gnz)
          real(RP) :: m
          integer :: i, n

          n = D%gnz

  ! initialize c-prime
          m = b1 + lambdas
          x(1) = x(1) / m
          if (n==1) return
          cp(1) = c1 / m
  ! solve for vectors c-prime and x-prime
          do i = 2, n-1
            m = D%mat_b(i) + lambdas - cp(i-1)*D%mat_a(i)
            cp(i) = D%mat_c(i)/m
            x(i) = (x(i)-x(i-1)*D%mat_a(i)) / m
          enddo
          x(n) = (x(n)-x(n-1)*D%mat_a(n)) / &
                 (D%mat_b(n) + lambdas - cp(n-1)*D%mat_a(n))
  ! solve for x from the vectors c-prime and x-prime
          do i = n-1, 1, -1
            x(i) = x(i) - cp(i)*x(i+1)
          end do

      end subroutine solve_tridiag 

    end subroutine solve_tridiagonal_1d_real_z
#endif
