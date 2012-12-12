program testpoisson
  use precisions
  use poisfft
  use vtkarray
  
  implicit none
  
  real(RP), parameter :: pi = 3.141592653589793238462_RP!pi = 4*atan(1._RP)!
  integer :: nx = 16, ny = 16, nz = 16
  real(RP), dimension(:,:,:), allocatable :: Phi, RHS
  real(RP), dimension(:,:), allocatable :: Phi2D, RHS2D
  real(RP), dimension(:), allocatable :: Phi1D, RHS1D
  real(RP) :: dx,dy,dz
  integer i,j,k,niters,inunit
  real(RP) :: R,x,y,z,tmp1,tmp2
  integer(8) :: t1,t2,trate
  type(PoisFFT_Solver3D) :: Solver
  type(PoisFFT_Solver2D) :: Solver2D
  type(PoisFFT_Solver1D) :: Solver1D
  character(len=5) :: ch5

  call system_clock(count_rate=trate)
  
!   read (*,*) nx,ny,nz
  if (command_argument_count()>=3) then
    call get_command_argument(1,value=ch5)
    read(ch5,'(i5)') nx
    call get_command_argument(2,value=ch5)
    read(ch5,'(i5)') ny
    call get_command_argument(3,value=ch5)
    read(ch5,'(i5)') nz
  endif
  
  if (command_argument_count()>=4) then
    call get_command_argument(4,value=ch5)
    read(ch5,'(i5)') niters
  else
    niters=10000
  endif
 
  dx=2*pi/nx;dy=2*pi/ny;dz=2*pi/nz
  
  allocate(RHS(nx,ny,nz))
  allocate(Phi(0:nx+1,0:ny+1,0:nz+1))

  allocate(RHS2D(nx,ny))
  allocate(Phi2D(0:nx+1,0:ny+1))

  allocate(RHS1D(nx))
  allocate(Phi1D(0:nx+1))


!   open(newunit=inunit,file="input.vtk")
!   read(unit=inunit,fmt="(13/)")
  do k=1,nz
   do j=1,ny
    do i=1,nx
     x=(i-1._RP/2)*dx
     y=(j-1._RP/2)*dy
     z=(k-1._RP/2)*dz
     !RHS(i,j,k)=cos(2*x)*cos(2*y)*cos(2*z)!(pi/2-abs(x-pi))*(pi/2-abs(y-pi))*(pi/2-abs(z-pi))!
     call random_number(RHS(i,j,k))
!      read(unit=inunit,fmt=*) RHS(i,j,k)
     call random_number(Phi(i,j,k))
    enddo
   enddo
  enddo
!   close(inunit)
  RHS=RHS-sum(RHS)/size(RHS)




  do j=1,ny
    do i=1,nx
     x=(i-1._RP/2)*dx
     y=(j-1._RP/2)*dy
     z=(k-1._RP/2)*dz
     !RHS(i,j,k)=cos(2*x)*cos(2*y)*cos(2*z)!(pi/2-abs(x-pi))*(pi/2-abs(y-pi))*(pi/2-abs(z-pi))!
     call random_number(RHS2D(i,j))
     call random_number(Phi2D(i,j))
    enddo
  enddo
  RHS2D=RHS2D-sum(RHS2D)/size(RHS2D)



  do i=1,nx
    x=(i-1._RP/2)*dx
!     RHS1D(i)=(pi/2-abs(x-pi))
    call random_number(RHS1D(i))
    call random_number(Phi1D(i))
  enddo
  RHS1D=RHS1D-sum(RHS1D)/size(RHS1D)
















  


  write (*,*) "3D PPNs"
  Phi(0,:,:)=Phi(nx,:,:)
  Phi(:,0,:)=Phi(:,ny,:)
  Phi(:,:,0)=Phi(:,:,1)
  Phi(nx+1,:,:)=Phi(1,:,:)
  Phi(:,ny+1,:)=Phi(:,1,:)
  Phi(:,:,nz+1)=Phi(:,:,nz)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,2,2],R)

  write (*,*) R

  Solver = PoisFFT_Solver3D_New(nx,ny,nz,dx,dy,dz,[(PoisFFT_PERIODIC, i=1,4),(PoisFFT_NeumannStag, i=5,6)])

  do i=1,3
    !call cpu_time(t1)
    call system_clock(count=t1)

    call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

    !call cpu_time(t2)
    call system_clock(count=t2)
    write(*,*) "solver cpu time", real(t2-t1)/real(trate)
  end do

  call PoisFFT_Solver3D_DeallocateData(Solver)

  Phi(0,:,:)=Phi(nx,:,:)
  Phi(:,0,:)=Phi(:,ny,:)
  Phi(:,:,0)=Phi(:,:,1)
  Phi(nx+1,:,:)=Phi(1,:,:)
  Phi(:,ny+1,:)=Phi(:,1,:)
  Phi(:,:,nz+1)=Phi(:,:,nz)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,2,2],R)
! 
  write (*,*) R


!  write(*,*) "sum", sum(Phi(1:nx,1:ny,1:nz))
! 
!   Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! 
! 
!  call VtkArraySimple("out.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))
! 
!   do i=1,niters/100
!     call GS3D(nx,ny,nz,100,Phi,RHS,&
!                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                 [3,3,3,3,2,2])
!   enddo
! 
!   call Res3D(nx,ny,nz,Phi,RHS,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                [3,3,3,3,2,2],R)
! 
!   write (*,*) R
! 
!  write(*,*) "sum", sum(Phi(1:nx,1:ny,1:nz))
! 
!  Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! 
! 
! call VtkArraySimple("gs.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))

!   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))

  write (*,*) "---------"


















  write(*,*) "3D staggered Dirichlet"

  do k=1,nz
   do j=1,ny
    do i=1,nx
     call random_number(Phi(i,j,k))
    enddo
   enddo
  enddo


!  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))


  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [1,1,1,1,1,1],R)

  write (*,*) R

  Solver = PoisFFT_Solver3D_New(nx,ny,nz,dx,dy,dz,[(PoisFFT_DirichletStag, i=1,6)])

  do i=1,3
    !call cpu_time(t1)
    call system_clock(count=t1)

    call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

    !call cpu_time(t2)
    call system_clock(count=t2)
    write(*,*) "solver cpu time", real(t2-t1)/real(trate)
  end do

  call PoisFFT_Solver3D_DeallocateData(Solver)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [1,1,1,1,1,1],R)

  write (*,*) R




!   call VtkArraySimple("outdir.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))
! 
! !   do k=1,nz
! !    do j=1,ny
! !     do i=1,nx
! !      call random_number(Phi(i,j,k))
! !     enddo
! !    enddo
! !   enddo
! 
!   do i=1,niters/100
!     call GS3D(nx,ny,nz,100,Phi,RHS,&
!                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                 [1,1,1,1,1,1])
!   call Res3D(nx,ny,nz,Phi,RHS,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                [1,1,1,1,1,1],R)
! 
!   write (*,*) R
!   enddo
! 
!   call Res3D(nx,ny,nz,Phi,RHS,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                [1,1,1,1,1,1],R)
! 
!   write (*,*) R
! 
!  call VtkArraySimple("gsdir.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))
! 
! !   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))
! 
 write(*,*) "--------"
















  write(*,*) "3D staggered Neumann:"

  do k=1,nz
   do j=1,ny
    do i=1,nx
     call random_number(Phi(i,j,k))
    enddo
   enddo
  enddo


!  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))


  Phi(0,:,:)=Phi(1,:,:)
  Phi(:,0,:)=Phi(:,1,:)
  Phi(:,:,0)=Phi(:,:,1)
  Phi(nx+1,:,:)=Phi(nx,:,:)
  Phi(:,ny+1,:)=Phi(:,ny,:)
  Phi(:,:,nz+1)=Phi(:,:,nz)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [2,2,2,2,2,2],R)

  write (*,*) R

  Solver = PoisFFT_Solver3D_New(nx,ny,nz,dx,dy,dz,[(PoisFFT_NeumannStag, i=1,6)])

  do i=1,3
    !call cpu_time(t1)
    call system_clock(count=t1)

    call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

    !call cpu_time(t2)
    call system_clock(count=t2)
    write(*,*) "solver cpu time", real(t2-t1)/real(trate)
  end do

  call PoisFFT_Solver3D_DeallocateData(Solver)

  Phi(0,:,:)=Phi(1,:,:)
  Phi(:,0,:)=Phi(:,1,:)
  Phi(:,:,0)=Phi(:,:,1)
  Phi(nx+1,:,:)=Phi(nx,:,:)
  Phi(:,ny+1,:)=Phi(:,ny,:)
  Phi(:,:,nz+1)=Phi(:,:,nz)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [2,2,2,2,2,2],R)

  write (*,*) R



! 
!   Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! 
!   tmp1=Phi(nx/2,ny/2,nz/2)
! 
!   call VtkArraySimple("outneum.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))
! 
! !   do k=1,nz
! !    do j=1,ny
! !     do i=1,nx
! !      call random_number(Phi(i,j,k))
! !     enddo
! !    enddo
! !   enddo
! 
!   do i=1,niters/100
!     call GS3D(nx,ny,nz,100,Phi,RHS,&
!                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                 [2,2,2,2,2,2])
!   enddo
! 
!   call Res3D(nx,ny,nz,Phi,RHS,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                [2,2,2,2,2,2],R)
! 
!   write (*,*) R
! 
!   Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! 
!   tmp2=Phi(nx/2,ny/2,nz/2)
! 
!  call VtkArraySimple("gsneum.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))

!   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))

 write(*,*) "--------"















  write(*,*) "3D Periodic"

  do k=1,nz
   do j=1,ny
    do i=1,nx
     call random_number(Phi(i,j,k))
    enddo
   enddo
  enddo



  Phi(0,:,:)=Phi(nx,:,:)
  Phi(:,0,:)=Phi(:,ny,:)
  Phi(:,:,0)=Phi(:,:,nz)
  Phi(nx+1,:,:)=Phi(1,:,:)
  Phi(:,ny+1,:)=Phi(:,1,:)
  Phi(:,:,nz+1)=Phi(:,:,1)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,3,3],R)

  write (*,*) R
  Solver = PoisFFT_Solver3D_New(nx,ny,nz,dx,dy,dz,[(PoisFFT_PERIODIC, i=1,6)])

  do i=1,3
    !call cpu_time(t1)
    call system_clock(count=t1)

    call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

    !call cpu_time(t2)
    call system_clock(count=t2)
    write(*,*) "solver cpu time", real(t2-t1)/real(trate)
  end do

  call PoisFFT_Solver3D_DeallocateData(Solver)

  Phi(0,:,:)=Phi(nx,:,:)
  Phi(:,0,:)=Phi(:,ny,:)
  Phi(:,:,0)=Phi(:,:,nz)
  Phi(nx+1,:,:)=Phi(1,:,:)
  Phi(:,ny+1,:)=Phi(:,1,:)
  Phi(:,:,nz+1)=Phi(:,:,1)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,3,3],R)

  write (*,*) R
!   Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! !
! call VtkArraySimple("periodic.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))
! 
!   do i=1,niters/100
!     call GS3D(nx,ny,nz,100,Phi,RHS,&
!                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                 [3,3,3,3,3,3])
!   enddo
! 
!   call Res3D(nx,ny,nz,Phi,RHS,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
!                [3,3,3,3,3,3],R)
! 
!   write (*,*) R
! 
!   Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! 
! 
! call VtkArraySimple("gsperiodic.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))
! 
!   call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))


 write(*,*) "--------"


















!   do j=1,ny
!     do i=1,nx
!       call random_number(Phi2D(i,j))
!     enddo
!   enddo
! 
! 
! !  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))
! 
! 
!   Phi2D(0,:)=Phi2D(1,:)
!   Phi2D(:,0)=Phi2D(:,1)
!   Phi2D(:,ny+1)=Phi2D(:,ny)
!   Phi2D(nx+1,:)=Phi2D(nx,:)
! 
! 
!   call Res2D(nx,ny,Phi2D,RHS2D,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),&
!                [2,2,2,2],R)
! 
!   write (*,*) R
! 
!   Solver2D = PoisFFT_Solver2D_New(nx,ny,dx,dy,[(PoisFFT_NeumannStag, i=1,4)])
! 
!   do i=1,3
!     !call cpu_time(t1)
!     call system_clock(count=t1)
! 
!     call Execute(Solver2D,Phi2D,RHS2D)
! 
!     !call cpu_time(t2)
!     call system_clock(count=t2)
!     write(*,*) "solver cpu time", real(t2-t1)/real(trate)
!   end do
! 
!   call DeallocateData(Solver2D)
! 
!   Phi2D(0,:)=Phi2D(1,:)
!   Phi2D(:,0)=Phi2D(:,1)
!   Phi2D(:,ny+1)=Phi2D(:,ny)
!   Phi2D(nx+1,:)=Phi2D(nx,:)
! 
!   call Res2D(nx,ny,Phi2D,RHS2D,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),&
!                [2,2,2,2],R)
! 
!   write (*,*) R
! 
! 
! 
! 
!   Phi2D=Phi2D-sum(Phi2D(1:nx,1:ny))/(nx*ny)
! 
!   tmp1=Phi2D(nx/2,ny/2)
! 
!   call VtkArraySimple("outneum2d.vtk",reshape(real(Phi2D(1:nx,1:ny),kind(1.)),[nx,ny,1]))
! 
! !   do k=1,nz
! !    do j=1,ny
! !     do i=1,nx
! !      call random_number(Phi(i,j,k))
! !     enddo
! !    enddo
! !   enddo
! 
!   do i=1,niters/100
!     call GS2D(nx,ny,100,Phi2D,RHS2D,&
!                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),&
!                 [2,2,2,2])
!   enddo
! 
!   call Res2D(nx,ny,Phi2D,RHS2D,&
!                dx**(-2),dx**(-2),dy**(-2),dy**(-2),&
!                [2,2,2,2],R)
! 
!   write (*,*) R
! 
!   Phi2D=Phi2D-sum(Phi2D(1:nx,1:ny))/(nx*ny)
! 
!   tmp2=Phi2D(nx/2,ny/2)
! 
!  call VtkArraySimple("gsneum2d.vtk",reshape(real(Phi2D(1:nx,1:ny),kind(1.)),[nx,ny,1]))
! 
! !   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))
! 
!  write(*,*) "--------"






















!   do i=1,nx
!      call random_number(Phi1D(i))
!   enddo
! 
! 
! !  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))
! 
! 
!   Phi1D(0)=Phi1D(1)
!   Phi1D(nx+1)=Phi1D(nx)
! 
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [3,3],R)
! 
!   write (*,*) R
! 
!   Solver1D = PoisFFT_Solver1D_New(nx,dx,[(PoisFFT_Periodic, i=1,2)])
! 
!   do i=1,3
!     !call cpu_time(t1)
!     call system_clock(count=t1)
! 
!     call Execute(Solver1D,Phi1D,RHS1D)
! 
!     !call cpu_time(t2)
!     call system_clock(count=t2)
!     write(*,*) "solver cpu time", real(t2-t1)/real(trate)
!   end do
! 
!   call DeallocateData(Solver1D)
! 
!   Phi1D(0)=Phi1D(1)
!   Phi1D(nx+1)=Phi1D(nx)
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [3,3],R)
! 
!   write (*,*) R
! 
! 
! 
! 
!   Phi1D=Phi1D-sum(Phi1D(1:nx))/(nx)
! 
!   tmp1=Phi1D(nx/2)
! 
!   call VtkArraySimple("outperiodic1d.vtk",reshape(real(Phi1D(1:nx),kind(1.)),[nx,1,1]))
! 
! !   do k=1,nz
! !    do j=1,ny
! !     do i=1,nx
! !      call random_number(Phi(i,j,k))
! !     enddo
! !    enddo
! !   enddo
! 
!   do i=1,niters/100
!     call GS1D(nx,100,Phi1D,RHS1D,&
!                 dx**(-2),dx**(-2),&
!                 [3,3])
!   enddo
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [3,3],R)
! 
!   write (*,*) R
! 
!   Phi1D=Phi1D-sum(Phi1D(1:nx))/(nx)
! 
!   tmp2=Phi1D(nx/2)
! 
!  call VtkArraySimple("gsperiodic1d.vtk",reshape(real(Phi1D(1:nx),kind(1.)),[nx,1,1]))
! 
! !   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))
! 
!  write(*,*) "--------"





















!   do i=1,nx
!      call random_number(Phi1D(i))
!   enddo
! 
! 
! !  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))
! 
! 
!   Phi1D(0)=Phi1D(1)
!   Phi1D(nx+1)=Phi1D(nx)
! 
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [1,1],R)
! 
!   write (*,*) R
! 
!   Solver1D = PoisFFT_Solver1D_New(nx,dx,[(PoisFFT_DirichletStag, i=1,2)])
! 
!   do i=1,3
!     !call cpu_time(t1)
!     call system_clock(count=t1)
! 
!     call Execute(Solver1D,Phi1D,RHS1D)
! 
!     !call cpu_time(t2)
!     call system_clock(count=t2)
!     write(*,*) "solver cpu time", real(t2-t1)/real(trate)
!   end do
! 
!   call DeallocateData(Solver1D)
! 
!   Phi1D(0)=Phi1D(1)
!   Phi1D(nx+1)=Phi1D(nx)
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [1,1],R)
! 
!   write (*,*) R
! 
! 
!   call VtkArraySimple("outdir1d.vtk",reshape(real(Phi1D(1:nx),kind(1.)),[nx,1,1]))
! 
! !   do k=1,nz
! !    do j=1,ny
! !     do i=1,nx
! !      call random_number(Phi(i,j,k))
! !     enddo
! !    enddo
! !   enddo
! 
!   do i=1,niters/100
!     call GS1D(nx,100,Phi1D,RHS1D,&
!                 dx**(-2),dx**(-2),&
!                 [1,1])
!   enddo
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [1,1],R)
! 
!   write (*,*) R
! 
! 
!  call VtkArraySimple("gsdir1d.vtk",reshape(real(Phi1D(1:nx),kind(1.)),[nx,1,1]))
! 
! !   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))
! 
!  write(*,*) "--------"












!   do i=1,nx
!      call random_number(Phi1D(i))
!   enddo
! 
! 
! !  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))
! 
! 
!   Phi1D(0)=Phi1D(1)
!   Phi1D(nx+1)=Phi1D(nx)
! 
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [2,2],R)
! 
!   write (*,*) R
! 
!   Solver1D = PoisFFT_Solver1D_New(nx,dx,[(PoisFFT_NeumannStag, i=1,2)])
! 
!   do i=1,3
!     !call cpu_time(t1)
!     call system_clock(count=t1)
! 
!     call Execute(Solver1D,Phi1D,RHS1D)
! 
!     !call cpu_time(t2)
!     call system_clock(count=t2)
!     write(*,*) "solver cpu time", real(t2-t1)/real(trate)
!   end do
! 
!   call DeallocateData(Solver1D)
! 
!   Phi1D(0)=Phi1D(1)
!   Phi1D(nx+1)=Phi1D(nx)
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [2,2],R)
! 
!   write (*,*) R
! 
! 
! 
! 
!   Phi1D=Phi1D-sum(Phi1D(1:nx))/(nx)
! 
!   tmp1=Phi1D(nx/2)
! 
!   call VtkArraySimple("outneum1d.vtk",reshape(real(Phi1D(1:nx),kind(1.)),[nx,1,1]))
! 
! !   do k=1,nz
! !    do j=1,ny
! !     do i=1,nx
! !      call random_number(Phi(i,j,k))
! !     enddo
! !    enddo
! !   enddo
! 
!   do i=1,niters/100
!     call GS1D(nx,100,Phi1D,RHS1D,&
!                 dx**(-2),dx**(-2),&
!                 [2,2])
!   enddo
! 
!   call Res1D(nx,Phi1D,RHS1D,&
!                dx**(-2),dx**(-2),&
!                [2,2],R)
! 
!   write (*,*) R
! 
!   Phi1D=Phi1D-sum(Phi1D(1:nx))/(nx)
! 
!   tmp2=Phi1D(nx/2)
! 
!  call VtkArraySimple("gsneum1d.vtk",reshape(real(Phi1D(1:nx),kind(1.)),[nx,1,1]))
! 
! !   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))
! 
!  write(*,*) "--------"









  write(*,*) "average of RHS:", sum(RHS1D)/(nx)

  contains
  
    subroutine Res1D(nx,Phi,RHS,&
                     Aw,Ae,&
                     Btype,R)

      implicit none
      intrinsic mod, abs, max

      integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
      integer, parameter :: Ea=1,We=2

      integer,intent(in)   :: nx
      real(KND),dimension(0:nx+1),intent(inout)::Phi
      real(KND),dimension(1:nx),intent(in)::RHS
      real(KND),intent(in) :: Aw,Ae
      integer,intent(in) :: Btype(2)
      real(KND),intent(out) :: R
      integer i,l
      real(KND) :: p,Ap

      R=0
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify(k,i), reduce(max:R)
           do i=1,nx
                p=0
                Ap=0
                if (i>1) then
                          p=p+Phi(i-1)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx)*Ae
                          Ap=Ap+Ae
                endif
  
                p=p-RHS(i)
                p=abs(-p +Ap*Phi(i))
                R=max(R,abs(p))

!                 write(*,'(x,g10.5)',advance="no") p
           enddo

!         write(*,*)
   endsubroutine Res1D










    subroutine Res2D(nx,ny,Phi,RHS,&
                     Aw,Ae,As,An,&
                     Btype,R)

      implicit none
      intrinsic mod, abs, max

      integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
      integer, parameter :: Ea=1,We=2,So=3,No=4

      integer,intent(in)   :: nx,ny
      real(KND),dimension(0:nx+1,0:ny+1),intent(inout)::Phi
      real(KND),dimension(1:nx,1:ny),intent(in)::RHS
      real(KND),intent(in) :: Aw,Ae
      real(KND),intent(in) :: As,An
      integer,intent(in) :: Btype(4)
      real(KND),intent(out) :: R
      integer i,j,k,l
      real(KND) :: p,Ap

      R=0
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify(k,i), reduce(max:R)
           do j=1,ny
              do i=1,nx
                p=0
                Ap=0
                if (i>1) then
                          p=p+Phi(i-1,j)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx,j)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1,j)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx,j)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny)*As
                          Ap=Ap+As
                 elseif (Btype(So)==DIRICHLET) then
                          p=p-Phi(i,1)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1)*An
                          Ap=Ap+An
                elseif (Btype(No)==DIRICHLET) then
                          p=p-Phi(i,ny)*An
                          Ap=Ap+An
                endif

                p=p-RHS(i,j)
                p=abs(-p +Ap*Phi(i,j))
                R=max(R,abs(p))

!                 write(*,'(x,g10.5)',advance="no") p
               enddo
           enddo
!         write(*,*)
   endsubroutine Res2D









    subroutine Res3D(nx,ny,nz,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype,R)

      implicit none
      intrinsic mod, abs, max

      integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
      integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6

      integer,intent(in)   :: nx,ny,nz
      real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
      real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
      real(KND),intent(in) :: Aw,Ae
      real(KND),intent(in) :: As,An
      real(KND),intent(in) :: Ab,At
      integer,intent(in) :: Btype(6)
      real(KND),intent(out) :: R
      integer i,j,k,l
      real(KND) :: p,Ap

      R=0
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify(k,i), reduce(max:R)
        do k=1,nz
          do j=1,ny
             do i=1,nx
                p=0
                Ap=0
                if (i>1) then
                          p=p+Phi(i-1,j,k)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx,j,k)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1,j,k)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j,k)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j,k)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx,j,k)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1,k)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny,k)*As
                          Ap=Ap+As
                 elseif (Btype(So)==DIRICHLET) then
                          p=p-Phi(i,1,k)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1,k)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1,k)*An
                          Ap=Ap+An
                elseif (Btype(No)==DIRICHLET) then
                          p=p-Phi(i,ny,k)*An
                          Ap=Ap+An
                endif
                if (k>1) then
                          p=p+Phi(i,j,k-1)*Ab
                          Ap=Ap+Ab
                elseif (Btype(Bo)==PERIODIC) then
                          p=p+Phi(i,j,nz)*Ab
                          Ap=Ap+Ab
                 elseif (Btype(Bo)==DIRICHLET) then
                          p=p-Phi(i,j,1)*Ab
                          Ap=Ap+Ab
                endif
                if (k<nz) then
                          p=p+Phi(i,j,k+1)*At
                          Ap=Ap+At
                elseif (Btype(To)==PERIODIC) then
                          p=p+Phi(i,j,1)*At
                          Ap=Ap+At
                elseif (Btype(To)==DIRICHLET) then
                          p=p-Phi(i,j,nz)*At
                          Ap=Ap+At
                endif
                p=p-RHS(i,j,k)
                p=abs(-p +Ap*Phi(i,j,k))
                R=max(R,abs(p))

!                 write(*,'(x,g10.5)',advance="no") p
               enddo
           enddo
        enddo
!         write(*,*)
   endsubroutine Res3D











  subroutine GS1D(nx,nit,Phi,RHS,&
                     Aw,Ae,&
                     Btype)
    implicit none

    integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
    integer, parameter :: Ea=1,We=2


    integer,intent(in)   :: nx,nit
    real(KND),dimension(0:nx+1),intent(inout)::Phi
    real(KND),dimension(1:nx),intent(in)::RHS
    real(KND),intent(in) :: Aw,Ae
    integer,intent(in) :: Btype(2)
    integer i,j,l
    real(KND) :: p,Ap
    intrinsic mod

    do l=1,nit
             do i=1,nx
               p=0
               Ap=0
                if (i>1) then
                          p=p+Phi(i-1)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx)*Ae
                          Ap=Ap+Ae
                endif

               p=p-RHS(i)

               p=p/Ap
               Phi(i)=Phi(i) + (p-Phi(i))
             enddo
    enddo

  endsubroutine GS1D





  subroutine GS2D(nx,ny,nit,Phi,RHS,&
                     Aw,Ae,As,An,&
                     Btype)
    implicit none

    integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
    integer, parameter :: Ea=1,We=2,So=3,No=4


    integer,intent(in)   :: nx,ny,nit
    real(KND),dimension(0:nx+1,0:ny+1),intent(inout)::Phi
    real(KND),dimension(1:nx,1:ny),intent(in)::RHS
    real(KND),intent(in) :: Aw,Ae
    real(KND),intent(in) :: As,An
    integer,intent(in) :: Btype(4)
    integer i,j,l
    real(KND) :: p,Ap
    intrinsic mod

    do l=1,nit
          do j=1,ny
              do i=1+mod(j+1,2),nx,2
               p=0
               Ap=0
                if (i>1) then
                          p=p+Phi(i-1,j)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx,j)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1,j)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx,j)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny)*As
                          Ap=Ap+As
                 elseif (Btype(So)==DIRICHLET) then
                          p=p-Phi(i,1)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1)*An
                          Ap=Ap+An
                elseif (Btype(No)==DIRICHLET) then
                          p=p-Phi(i,ny)*An
                          Ap=Ap+An
                endif

               p=p-RHS(i,j)

               p=p/Ap
               Phi(i,j)=p
              enddo
          enddo
          do j=1,ny
              do i=1+mod(j,2),nx,2
               p=0
               Ap=0
                if (i>1) then
                          p=p+Phi(i-1,j)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx,j)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1,j)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx,j)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny)*As
                          Ap=Ap+As
                 elseif (Btype(So)==DIRICHLET) then
                          p=p-Phi(i,1)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1)*An
                          Ap=Ap+An
                elseif (Btype(No)==DIRICHLET) then
                          p=p-Phi(i,ny)*An
                          Ap=Ap+An
                endif

               p=p-RHS(i,j)

               p=p/Ap
               Phi(i,j)=p
              enddo
          enddo
    enddo

  endsubroutine GS2D





  subroutine GS3D(nx,ny,nz,nit,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype)
     implicit none

     integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
     integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6


     integer,intent(in)   :: nx,ny,nz,nit
     real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
     real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
     real(KND),intent(in) :: Aw,Ae
     real(KND),intent(in) :: As,An
     real(KND),intent(in) :: Ab,At
     integer,intent(in) :: Btype(6)
     integer i,j,k,l
     real(KND) :: p,Ap
     intrinsic mod

    do l=1,nit
       do k=1,nz
          do j=1,ny
              do i=1,nx
               p=0
               Ap=0
                if (i>1) then
                          p=p+Phi(i-1,j,k)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx,j,k)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==DIRICHLET) then
                          p=p-Phi(1,j,k)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j,k)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j,k)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==DIRICHLET) then
                          p=p-Phi(nx,j,k)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1,k)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny,k)*As
                          Ap=Ap+As
                 elseif (Btype(So)==DIRICHLET) then
                          p=p-Phi(i,1,k)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1,k)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1,k)*An
                          Ap=Ap+An
                elseif (Btype(No)==DIRICHLET) then
                          p=p-Phi(i,ny,k)*An
                          Ap=Ap+An
                endif
                if (k>1) then
                          p=p+Phi(i,j,k-1)*Ab
                          Ap=Ap+Ab
                elseif (Btype(Bo)==PERIODIC) then
                          p=p+Phi(i,j,nz)*Ab
                          Ap=Ap+Ab
                 elseif (Btype(Bo)==DIRICHLET) then
                          p=p-Phi(i,j,1)*Ab
                          Ap=Ap+Ab
                endif
                if (k<nz) then
                          p=p+Phi(i,j,k+1)*At
                          Ap=Ap+At
                elseif (Btype(To)==PERIODIC) then
                          p=p+Phi(i,j,1)*At
                          Ap=Ap+At
                elseif (Btype(To)==DIRICHLET) then
                          p=p-Phi(i,j,nz)*At
                          Ap=Ap+At
                endif
               p=p-RHS(i,j,k)

               p=p/Ap
               Phi(i,j,k)=p
          enddo
        enddo
      enddo
    enddo

  endsubroutine GS3D

end program testpoisson
