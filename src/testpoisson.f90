program testpoisson

  use poisfft
  use vtkarray
  
  implicit none
  
  real(RP), parameter :: pi = 4*atan(1._RP)!3.141592653589793238462_RP
  integer :: nx = 16, ny = 16, nz = 16
  real(RP), dimension(:,:,:), allocatable :: Phi, RHS
  real(RP) :: dx,dy,dz

  integer i,j,k
  real(RP) :: R,x,y,z,t1,t2,tmp1,tmp2

  type(PoisFFT_Solver3D) :: Solver

  read (*,*) nx,ny,nz

  dx=2*pi/nx;dy=2*pi/ny;dz=2*pi/nz
  
  allocate(RHS(nx,ny,nz))
  allocate(Phi(0:nx+1,0:ny+1,0:nz+1))

  do k=1,nz
   do j=1,ny
    do i=1,nx
     x=(i-1._RP/2)*dx
     y=(j-1._RP/2)*dy
     z=(k-1._RP/2)*dz
     RHS(i,j,k)=cos(2*x)*cos(2*y)*cos(2*z)
     call random_number(Phi(i,j,k))
    enddo
   enddo
  enddo


  call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))
  

  Phi(0,:,:)=Phi(1,:,:)
  Phi(:,0,:)=Phi(:,1,:)
  Phi(:,:,0)=Phi(:,:,1)
  Phi(nx+1,:,:)=Phi(nx,:,:)
  Phi(:,ny+1,:)=Phi(:,ny,:)
  Phi(:,:,nz+1)=Phi(:,:,nz)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,1,1],R)

  write (*,*) R

  Solver = PoisFFT_Solver3D_New(nx,ny,nz,dx,dy,dz,[(PoisFFT_PERIODIC, i=1,4),(PoisFFT_NeumannStag, i=5,6)])

  call cpu_time(t1)

  call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

  call cpu_time(t2)
  write(*,*) "solver cpu time", t2-t1

  call cpu_time(t1)

  call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

  call cpu_time(t2)
  write(*,*) "solver cpu time", t2-t1

  call cpu_time(t1)

  call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

  call cpu_time(t2)
  write(*,*) "solver cpu time", t2-t1

  call PoisFFT_Solver3D_DeallocateData(Solver)

  Phi(0,:,:)=Phi(1,:,:)
  Phi(:,0,:)=Phi(:,1,:)
  Phi(:,:,0)=Phi(:,:,1)
  Phi(nx+1,:,:)=Phi(nx,:,:)
  Phi(:,ny+1,:)=Phi(:,ny,:)
  Phi(:,:,nz+1)=Phi(:,:,nz)

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,1,1],R)

  write (*,*) R




  Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
  
  tmp1=Phi(nx/2,ny/2,nz/2)
  
  call VtkArraySimple("out.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))

  do i=1,100
    call GS3D(nx,ny,nz,100,Phi,RHS,&
                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
                [3,3,3,3,1,1])
  enddo

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,1,1],R)

  write (*,*) R

  Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)

  tmp2=Phi(nx/2,ny/2,nz/2)

  call VtkArraySimple("gs.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))

!   call VtkArraySimple("rhs2.vtk",real(RHS,kind(1.)))

write (*,*) "ratio:", tmp1/tmp2





















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

  call cpu_time(t1)

  call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

  call cpu_time(t2)
  write(*,*) "solver cpu time", t2-t1

  call cpu_time(t1)

  call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

  call cpu_time(t2)
  write(*,*) "solver cpu time", t2-t1

  call cpu_time(t1)

  call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

  call cpu_time(t2)
  write(*,*) "solver cpu time", t2-t1

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
  Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)
! 
  call VtkArraySimple("periodic.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))

  do i=1,30
    call GS3D(nx,ny,nz,100,Phi,RHS,&
                dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
                [3,3,3,3,3,3])
  enddo

  call Res3D(nx,ny,nz,Phi,RHS,&
               dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
               [3,3,3,3,3,3],R)

  write (*,*) R

  Phi=Phi-sum(Phi(1:nx,1:ny,1:nz))/(nx*ny*nz)


  call VtkArraySimple("gsperiodic.vtk",real(Phi(1:nx,1:ny,1:nz),kind(1.)))

!   call VtkArraySimple("rhs.vtk",real(RHS,kind(1.)))

  contains
  
    subroutine Res1D(nx,Phi,RHS,&
                     Aw,Ae,&
                     Btype,R)

      implicit none
      intrinsic mod, abs, max

      integer,parameter:: KND=RP,PERIODIC=3
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
                endif
                if (i<nx) then
                          p=p+Phi(i+1)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1)*Ae
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

      integer,parameter:: KND=RP,PERIODIC=3
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
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1)*An
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

      integer,parameter:: KND=RP,PERIODIC=3
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
           do i=1,nx
               do j=1,ny
                p=0
                Ap=0
                if (i>1) then
                          p=p+Phi(i-1,j,k)*Aw
                          Ap=Ap+Aw
                elseif (Btype(We)==PERIODIC) then
                          p=p+Phi(nx,j,k)*Aw
                          Ap=Ap+Aw
                endif
                if (i<nx) then
                          p=p+Phi(i+1,j,k)*Ae
                          Ap=Ap+Ae
                elseif (Btype(Ea)==PERIODIC) then
                          p=p+Phi(1,j,k)*Ae
                          Ap=Ap+Ae
                endif
                if (j>1) then
                          p=p+Phi(i,j-1,k)*As
                          Ap=Ap+As
                elseif (Btype(So)==PERIODIC) then
                          p=p+Phi(i,ny,k)*As
                          Ap=Ap+As
                endif
                if (j<ny) then
                          p=p+Phi(i,j+1,k)*An
                          Ap=Ap+An
                elseif (Btype(No)==PERIODIC) then
                          p=p+Phi(i,1,k)*An
                          Ap=Ap+An
                endif
                if (k>1) then
                          p=p+Phi(i,j,k-1)*Ab
                          Ap=Ap+Ab
                elseif (Btype(Bo)==PERIODIC) then
                          p=p+Phi(i,j,nz)*Ab
                          Ap=Ap+Ab
                endif
                if (k<nz) then
                          p=p+Phi(i,j,k+1)*At
                          Ap=Ap+At
                elseif (Btype(To)==PERIODIC) then
                          p=p+Phi(i,j,1)*At
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

   subroutine GS3D(nx,ny,nz,nit,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype)
     implicit none

     integer,parameter:: KND=RP,PERIODIC=3
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
     !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i)
       do k=1,nz
          do i=1,nx
              do j=1+mod(i+k,2),ny,2
               p=0
               Ap=0
               if (i>1) then
                         p=p+Phi(i-1,j,k)*Aw
                         Ap=Ap+Aw
               elseif (Btype(We)==PERIODIC) then
                         p=p+Phi(nx,j,k)*Aw
                         Ap=Ap+Aw
               endif
               if (i<nx) then
                         p=p+Phi(i+1,j,k)*Ae
                         Ap=Ap+Ae
               elseif (Btype(Ea)==PERIODIC) then
                         p=p+Phi(1,j,k)*Ae
                         Ap=Ap+Ae
               endif
               if (j>1) then
                         p=p+Phi(i,j-1,k)*As
                         Ap=Ap+As
               elseif (Btype(So)==PERIODIC) then
                         p=p+Phi(i,ny,k)*As
                         Ap=Ap+As
               endif
               if (j<ny) then
                         p=p+Phi(i,j+1,k)*An
                         Ap=Ap+An
               elseif (Btype(No)==PERIODIC) then
                         p=p+Phi(i,1,k)*An
                         Ap=Ap+An
               endif
               if (k>1) then
                         p=p+Phi(i,j,k-1)*Ab
                         Ap=Ap+Ab
               elseif (Btype(Bo)==PERIODIC) then
                         p=p+Phi(i,j,nz)*Ab
                         Ap=Ap+Ab
               endif
               if (k<nz) then
                         p=p+Phi(i,j,k+1)*At
                         Ap=Ap+At
               elseif (Btype(To)==PERIODIC) then
                         p=p+Phi(i,j,1)*At
                         Ap=Ap+At
               endif
               p=p-RHS(i,j,k)

               p=p/Ap
               Phi(i,j,k)=p
              enddo
          enddo
      enddo
    !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i)
       do k=1,nz
          do i=1,nx
              do j=1+mod(i+k+1,2),ny,2
               p=0
               Ap=0
               if (i>1) then
                         p=p+Phi(i-1,j,k)*Aw
                         Ap=Ap+Aw
               elseif (Btype(We)==PERIODIC) then
                         p=p+Phi(nx,j,k)*Aw
                         Ap=Ap+Aw
               endif
               if (i<nx) then
                         p=p+Phi(i+1,j,k)*Ae
                         Ap=Ap+Ae
               elseif (Btype(Ea)==PERIODIC) then
                         p=p+Phi(1,j,k)*Ae
                         Ap=Ap+Ae
               endif
               if (j>1) then
                         p=p+Phi(i,j-1,k)*As
                         Ap=Ap+As
               elseif (Btype(So)==PERIODIC) then
                         p=p+Phi(i,ny,k)*As
                         Ap=Ap+As
               endif
               if (j<ny) then
                         p=p+Phi(i,j+1,k)*An
                         Ap=Ap+An
               elseif (Btype(No)==PERIODIC) then
                         p=p+Phi(i,1,k)*An
                         Ap=Ap+An
               endif
               if (k>1) then
                         p=p+Phi(i,j,k-1)*Ab
                         Ap=Ap+Ab
               elseif (Btype(Bo)==PERIODIC) then
                         p=p+Phi(i,j,nz)*Ab
                         Ap=Ap+Ab
               endif
               if (k<nz) then
                         p=p+Phi(i,j,k+1)*At
                         Ap=Ap+At
               elseif (Btype(To)==PERIODIC) then
                         p=p+Phi(i,j,1)*At
                         Ap=Ap+At
               endif
               p=p-RHS(i,j,k)

               p=p/Ap
               Phi(i,j,k)=p
              enddo
          enddo
      enddo
    enddo
      Phi=Phi-Phi(nx/2,ny/2,nz/2)
  endsubroutine GS3D

end program testpoisson
