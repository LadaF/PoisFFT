      !$omp do
      do k=2,D%nz
        do j=2,D%ny
          do i=2,D%nx
            D%xwork(i,j,k) = D%xwork(i,j,k) / (denomx(i) + denomy(j) + denomz(k))
          end do
        end do
      end do
      !$omp end do
      
      !$omp do
        do j=2,D%ny
          do i=2,D%nx
            D%xwork(i,j,1) = D%xwork(i,j,1) / (denomx(i) + denomy(j))
          end do
        end do
      !$omp end do

      !$omp do
      do k=2,D%nz
          do i=2,D%nx
            D%xwork(i,1,k) = D%xwork(i,1,k) / (denomx(i) + denomz(k))
          end do
      end do
      !$omp end do

      !$omp do
      do k=2,D%nz
        do j=2,D%ny
            D%xwork(1,j,k) = D%xwork(1,j,k) / (denomy(j) + denomz(k))
        end do
      end do
      !$omp end do
      
      !$omp do
          do i=2,D%nx
            D%xwork(i,1,1) = D%xwork(i,1,1) / (denomx(i))
          end do
      !$omp end do

      !$omp do
        do j=2,D%ny
            D%xwork(1,j,1) = D%xwork(1,j,1) / (denomy(j))
        end do
      !$omp end do

      !$omp do
      do k=2,D%nz
            D%xwork(1,1,k) = D%xwork(1,1,k) / (denomz(k))
      end do
      !$omp end do

      