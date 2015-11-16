module Endianness

  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd

  logical,save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd32
    module procedure BigEnd64
  end interface

  contains

    subroutine GetEndianness
      integer(int8),dimension(4):: bytes !may not work on some processors

      bytes=transfer(1_int32,bytes,4)
      if (bytes(4)==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

    elemental function BigEnd32(x) result(res)
      real(real32) :: res
      real(real32),intent(in)::x
      integer(int8),dimension(4):: bytes !may not work on some processors

      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes,4)
        res = transfer(bytes(4:1:-1),res)
      endif
    end function BigEnd32

    elemental function BigEnd64(x) result(res)
      real(real64) :: res
      real(real64),intent(in)::x
      integer(int8),dimension(8):: bytes !may not work on some processors

      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes,8)
        res = transfer(bytes(8:1:-1),res)
      endif
    end function BigEnd64

end module Endianness

module parameters
   use PoisFFT_Precisions
   use iso_c_binding
   integer,parameter :: RP = DRP
   real(RP), parameter :: pi = 3.141592653589793238462_RP!pi = 4*atan(1._RP)!
end module

module my_mpi
  use parameters
  use mpi
  
  integer :: nims, npxyz(3), pxyz(3)
  integer :: nxims, nyims, nzims
  integer :: myim, myrank, iim, jim, kim
  integer :: w_im, e_im, s_im, n_im, b_im, t_im
  integer :: w_rank, e_rank, s_rank, n_rank, b_rank, t_rank
  integer :: glob_comm, cart_comm, cart_comm_dim = -1
  logical :: master = .false.
  integer :: MPI_RP = -huge(1)
  
  interface error_stop
    module procedure error_stop_int
    module procedure error_stop_char
  end interface
  
 contains

  integer function this_image() result(res)
    integer ie
    call MPI_Comm_rank(glob_comm, res, ie)
    res = res + 1
    if (ie/=0) call error_stop("MPI_Comm_rank ERROR")
  end function

  integer function num_images() result(res)
    integer ie
    call MPI_Comm_size(glob_comm, res, ie)
    if (ie/=0) call error_stop("MPI_Comm_size ERROR")
  end function

  integer function image_index(sub) result(res)
    integer, intent(in) :: sub(3)
    integer ie
    
    if (cart_comm_dim==-1) then
      call MPI_Cartdim_get(cart_comm, cart_comm_dim, ie)
      if (ie/=0) call error_stop("MPI_Cartdim_get")
    end if
    call MPI_Cart_rank(cart_comm, sub(3:4-cart_comm_dim:-1)-1, res, ie)  
    if (ie/=0) call error_stop("MPI_Cart_rank")
    res = res + 1
  end function
  
  subroutine get_image_coords()
    call MPI_Cart_coords(cart_comm, myrank, 3, pxyz, ie)
    if (ie/=0) call error_stop("MPI_Cart_coords")
        
     pxyz = pxyz(3:1:-1)
     
     iim = pxyz(1) + 1
     jim = pxyz(2) + 1
     kim = pxyz(3) + 1
     
     if (iim>1) then
       w_im = image_index([iim-1, jim, kim])
     else
       w_im = image_index([nxims, jim, kim])
     end if
     if (iim<nxims) then
       e_im = image_index([iim+1, jim, kim])
     else
       e_im = image_index([1, jim, kim])
     end if
     
     if (jim>1) then
       s_im = image_index([iim, jim-1, kim])
     else
       s_im = image_index([iim, nyims, kim])
     end if
     if (jim<nyims) then
       n_im = image_index([iim, jim+1, kim])
     else
       n_im = image_index([iim, 1, kim])
     end if

     if (kim>1) then
       b_im = image_index([iim, jim, kim-1])
     else
       b_im = image_index([iim, jim, nzims])
     end if
     if (kim<nzims) then
       t_im = image_index([iim, jim, kim+1])
     else
       t_im = image_index([iim, jim, 1])
     end if
     
     w_rank = w_im - 1
     e_rank = e_im - 1
     s_rank = s_im - 1
     n_rank = n_im - 1
     b_rank = b_im - 1
     t_rank = t_im - 1
  end subroutine
  
    
  subroutine error_stop_int(n)
    integer,intent(in) :: n
    integer ie
    if (master) then
      write(*,*) "ERROR",n
    end if
    call MPI_finalize(ie)
    stop
  end subroutine

  subroutine error_stop_char(ch)
    character(*),intent(in) :: ch
    integer ie
    if (master) then
      write(*,*) "ERROR",ch
    end if
    call MPI_finalize(ie)
    stop
  end subroutine

end module




module subs
   use parameters
   use my_mpi
   use PoisFFT
   
   implicit none
   

contains

    subroutine exchange_boundaries_3D(comm, Phi, nx, ny, nz, BCs)
      integer, intent(in) :: comm
      real(RP), intent(inout),contiguous :: Phi(0:,0:,0:)
      integer, intent(in) :: nx, ny, nz
      integer, intent(in) :: BCs(6)
      logical :: oddx, oddy, oddz, evenx, eveny, evenz
      integer ierr,tag,status(MPI_STATUS_SIZE)
      
      oddx = mod(iim,2) == 1
      evenx = .not. oddx
      
      oddy = mod(jim,2) == 1
      eveny = .not. oddy
      
      oddz = mod(kim,2) == 1
      evenz = .not. oddz
      

      call MPI_Barrier(comm, ierr)

      !internal boundaries
      if (oddx) then
        call send_w
      else
        call recv_e
      end if
      if (evenx) then
        call send_w
      else
        call recv_e
      end if

      if (oddx) then
        call send_e
      else
        call recv_w
      end if
      if (evenx) then
        call send_e
      else
        call recv_w
      end if

      if (oddy) then
        call send_s
      else
        call recv_n
      end if
      if (eveny) then
        call send_s
      else
        call recv_n
      end if

      if (oddy) then
        call send_n
      else
        call recv_s
      end if
      if (eveny) then
        call send_n
      else
        call recv_s
      end if

      if (oddz) then
        call send_b
      else
        call recv_t
      end if
      if (evenz) then
        call send_b
      else
        call recv_t
      end if

      if (oddz) then
        call send_t
      else
        call recv_b
      end if
      if (evenz) then
        call send_t
      else
        call recv_b
      end if
      

      !global domain boundaries
      if (BCs(1)==PoisFFT_PERIODIC) then
        if (nxims>1) then
          if (iim==1) then
            call send(Phi(1,1:ny,1:nz), w_rank)
          else if (iim==nxims) then
            call recv(Phi(nx+1,1:ny,1:nz), e_rank)
          end if
          if (iim==nxims) then
            call send(Phi(nx,1:ny,1:nz), e_rank)
          else if (iim==1) then
            call recv(Phi(0,1:ny,1:nz), w_rank)
          end if
        else
          Phi(0,:,:)=Phi(nx,:,:)
          Phi(nx+1,:,:)=Phi(1,:,:)
        end if
      else
        if (BCs(1)==PoisFFT_NeumannStag) then
          if (iim==1) Phi(0,:,:) = Phi(1,:,:)
        elseif (BCs(1)==PoisFFT_DirichletStag) then
          if (iim==1) Phi(0,:,:) = -Phi(1,:,:)
        end if
        if (BCs(2)==PoisFFT_NeumannStag) then
          if (iim==nxims) Phi(nx+1,:,:) = Phi(nx,:,:)
        elseif (BCs(2)==PoisFFT_DirichletStag) then
          if (iim==nxims) Phi(nx+1,:,:) = -Phi(nx,:,:)
        end if
      end if

        
      if (BCs(3)==PoisFFT_PERIODIC) then
        if (nyims>1) then
          if (jim==1) then
            call send(Phi(1:nx,1,1:nz), s_rank)
          else if (jim==nyims) then
            call recv(Phi(1:nx,ny+1,1:nz), n_rank)
          end if
          if (jim==nyims) then
            call send(Phi(1:nx,ny,1:nz), n_rank)
          else if (jim==1) then
            call recv(Phi(1:nx,0,1:nz), s_rank)
          end if
        else
          Phi(:,0,:)=Phi(:,ny,:)
          Phi(:,ny+1,:)=Phi(:,1,:)
        end if
      else
        if (BCs(3)==PoisFFT_NeumannStag) then
          if (jim==1) Phi(:,0,:) = Phi(:,1,:)
        elseif (BCs(3)==PoisFFT_DirichletStag) then
          if (jim==1) Phi(:,0,:) = -Phi(:,1,:)
        end if
        if (BCs(4)==PoisFFT_NeumannStag) then
          if (jim==nyims) Phi(:,ny+1,:) = Phi(:,ny,:)
        elseif (BCs(4)==PoisFFT_DirichletStag) then
          if (jim==nyims) Phi(:,ny+1,:) = -Phi(:,ny,:)
        end if
      end if

      if (BCs(5)==PoisFFT_PERIODIC) then
        if (nzims>1) then
          if (kim==1) then
            call send(Phi(1:nx,1:ny,1), b_rank)
          else if (kim==nzims) then
            call recv(Phi(1:nx,1:ny,nz+1), t_rank)
          end if
          if (kim==nzims) then
            call send(Phi(1:nx,1:ny,nz), t_rank)
          else if (kim==1) then
            call recv(Phi(1:nx,1:ny,0), b_rank)
          end if
        else
          Phi(:,:,0)=Phi(:,:,nz)
          Phi(:,:,nz+1)=Phi(:,:,1)
        end if
      else
        if (BCs(5)==PoisFFT_NeumannStag) then
          if (kim==1) Phi(:,:,0) = Phi(:,:,1)
        elseif (BCs(5)==PoisFFT_DirichletStag) then
          if (kim==1) Phi(:,:,0) = -Phi(:,:,1)
        end if
        if (BCs(6)==PoisFFT_NeumannStag) then
          if (kim==nzims) Phi(:,:,nz+1) = Phi(:,:,nz)
        elseif (BCs(6)==PoisFFT_DirichletStag) then
          if (kim==nzims) Phi(:,:,nz+1) = -Phi(:,:,nz)
        end if
      end if
              
      call MPI_Barrier(comm, ierr)
    
    
    contains
    
    
      subroutine send(a,to)
        real(RP), intent(in) :: a(:,:)
        integer, intent(in) :: to

        call MPI_Send(a, size(a) , MPI_RP, to, 1, comm, ierr)
        if (ierr/=0) stop "error sending MPI message."
      end subroutine
      
      subroutine recv(a,from)
        real(RP), intent(out) :: a(:,:)
        integer, intent(in) :: from

        call MPI_Recv(a, size(a) , MPI_RP, from, 1, comm, status, ierr)
        if (ierr/=0) stop "error sending MPI message."
      end subroutine
      

      subroutine send_w
        if (iim>1) then
          call send(Phi(1,1:ny,1:nz), w_rank)
        end if
      end subroutine
      subroutine recv_w
        if (iim>1) then
          call recv(Phi(0,1:ny,1:nz), w_rank)
        end if
      end subroutine
      subroutine send_e
        if (iim<nxims) then
          call send(Phi(nx,1:ny,1:nz), e_rank)
        end if
      end subroutine       
      subroutine recv_e
        if (iim<nxims) then
          call recv(Phi(nx+1,1:ny,1:nz), e_rank)
        end if
       end subroutine
      subroutine send_s
        if (jim>1) then
          call send(Phi(1:nx,1,1:nz), s_rank)
        end if
      end subroutine
      subroutine recv_s
        if (jim>1) then
          call recv(Phi(1:nx,0,1:nz), s_rank)
        end if
      end subroutine
      subroutine send_n
        if (jim<nyims) then
          call send(Phi(1:nx,ny,1:nz), n_rank)
        end if
      end subroutine
      subroutine recv_n
        if (jim<nyims) then
          call recv(Phi(1:nx,ny+1,1:nz), n_rank)
        end if
      end subroutine
      subroutine send_b
        if (kim>1) then
          call send(Phi(1:nx,1:ny,1), b_rank)
        end if
      end subroutine
      subroutine recv_b
        if (kim>1) then
          call recv(Phi(1:nx,1:ny,0), b_rank)
        end if
      end subroutine
      subroutine send_t
        if (kim<nzims) then
          call send(Phi(1:nx,1:ny,nz), t_rank)
        end if
      end subroutine
      subroutine recv_t
        if (kim<nzims) then
          call recv(Phi(1:nx,1:ny,nz+1), t_rank)
        end if
      end subroutine
      
    end subroutine exchange_boundaries_3D

  
  



    subroutine exchange_boundaries_2D(comm, Phi, nx, ny, BCs)
      integer, intent(in) :: comm
      real(RP), intent(inout),contiguous :: Phi(0:,0:)
      integer, intent(in) :: nx, ny
      integer, intent(in) :: BCs(4)
      logical :: oddx, oddy, evenx, eveny
      integer ierr,tag,status(MPI_STATUS_SIZE)
      
      oddx = mod(iim,2) == 1
      evenx = .not. oddx
      
      oddy = mod(jim,2) == 1
      eveny = .not. oddy
      

      call MPI_Barrier(comm, ierr)

      !internal boundaries
      if (oddx) then
        call send_w
      else
        call recv_e
      end if
      if (evenx) then
        call send_w
      else
        call recv_e
      end if

      if (oddx) then
        call send_e
      else
        call recv_w
      end if
      if (evenx) then
        call send_e
      else
        call recv_w
      end if

      if (oddy) then
        call send_s
      else
        call recv_n
      end if
      if (eveny) then
        call send_s
      else
        call recv_n
      end if

      if (oddy) then
        call send_n
      else
        call recv_s
      end if
      if (eveny) then
        call send_n
      else
        call recv_s
      end if


      !global domain boundaries
      if (BCs(1)==PoisFFT_PERIODIC) then
        if (nxims>1) then
          if (iim==1) then
            call send(Phi(1,1:ny), w_rank)
          else if (iim==nxims) then
            call recv(Phi(nx+1,1:ny), e_rank)
          end if
          if (iim==nxims) then
            call send(Phi(nx,1:ny), e_rank)
          else if (iim==1) then
            call recv(Phi(0,1:ny), w_rank)
          end if
        else
          Phi(0,:)=Phi(nx,:)
          Phi(nx+1,:)=Phi(1,:)
        end if
      else
        if (BCs(1)==PoisFFT_NeumannStag) then
          if (iim==1) Phi(0,:) = Phi(1,:)
        elseif (BCs(1)==PoisFFT_DirichletStag) then
          if (iim==1) Phi(0,:) = -Phi(1,:)
        end if
        if (BCs(2)==PoisFFT_NeumannStag) then
          if (iim==nxims) Phi(nx+1,:) = Phi(nx,:)
        elseif (BCs(2)==PoisFFT_DirichletStag) then
          if (iim==nxims) Phi(nx+1,:) = -Phi(nx,:)
        end if
      end if

        
      if (BCs(3)==PoisFFT_PERIODIC) then
        if (nyims>1) then
          if (jim==1) then
            call send(Phi(1:nx,1), s_rank)
          else if (jim==nyims) then
            call recv(Phi(1:nx,ny+1), n_rank)
          end if
          if (jim==nyims) then
            call send(Phi(1:nx,ny), n_rank)
          else if (jim==1) then
            call recv(Phi(1:nx,0), s_rank)
          end if
        else
          Phi(:,0)=Phi(:,ny)
          Phi(:,ny+1)=Phi(:,1)
        end if
      else
        if (BCs(3)==PoisFFT_NeumannStag) then
          if (jim==1) Phi(:,0) = Phi(:,1)
        elseif (BCs(3)==PoisFFT_DirichletStag) then
          if (jim==1) Phi(:,0) = -Phi(:,1)
        end if
        if (BCs(4)==PoisFFT_NeumannStag) then
          if (jim==nyims) Phi(:,ny+1) = Phi(:,ny)
        elseif (BCs(4)==PoisFFT_DirichletStag) then
          if (jim==nyims) Phi(:,ny+1) = -Phi(:,ny)
        end if
      end if

              
      call MPI_Barrier(comm, ierr)
    
    
    contains
    
    
      subroutine send(a,to)
        real(RP), intent(in) :: a(:)
        integer, intent(in) :: to

        call MPI_Send(a, size(a) , MPI_RP, to, 1, comm, ierr)
        if (ierr/=0) stop "error sending MPI message."
      end subroutine
      
      subroutine recv(a,from)
        real(RP), intent(out) :: a(:)
        integer, intent(in) :: from

        call MPI_Recv(a, size(a) , MPI_RP, from, 1, comm, status, ierr)
        if (ierr/=0) stop "error sending MPI message."
      end subroutine
      

      subroutine send_w
        if (iim>1) then
          call send(Phi(1,1:ny), w_rank)
        end if
      end subroutine
      subroutine recv_w
        if (iim>1) then
          call recv(Phi(0,1:ny), w_rank)
        end if
      end subroutine
      subroutine send_e
        if (iim<nxims) then
          call send(Phi(nx,1:ny), e_rank)
        end if
      end subroutine       
      subroutine recv_e
        if (iim<nxims) then
          call recv(Phi(nx+1,1:ny), e_rank)
        end if
       end subroutine
      subroutine send_s
        if (jim>1) then
          call send(Phi(1:nx,1), s_rank)
        end if
      end subroutine
      subroutine recv_s
        if (jim>1) then
          call recv(Phi(1:nx,0), s_rank)
        end if
      end subroutine
      subroutine send_n
        if (jim<nyims) then
          call send(Phi(1:nx,ny), n_rank)
        end if
      end subroutine
      subroutine recv_n
        if (jim<nyims) then
          call recv(Phi(1:nx,ny+1), n_rank)
        end if
      end subroutine

      
    end subroutine exchange_boundaries_2D

  
  




    subroutine Res3D(nx,ny,nz,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     R)

      implicit none
      intrinsic mod, abs, max

      integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
      integer, parameter :: We=1,Ea=2,So=3,No=4,Bo=5,To=6

      integer,intent(in)   :: nx,ny,nz
      real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
      real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
      real(KND),intent(in) :: Aw,Ae
      real(KND),intent(in) :: As,An
      real(KND),intent(in) :: Ab,At
      real(KND),intent(out) :: R
      integer i,j,k,l
      real(KND) :: p,Ap
      integer ie

      R=0
      Ap = Aw + Ae + As + An + Ab + At
      
      do k=1,nz
        do j=1,ny
           do i=1,nx
              p=0
              p=p+Phi(i-1,j,k)*Aw
              p=p+Phi(i+1,j,k)*Ae
              p=p+Phi(i,j-1,k)*As
              p=p+Phi(i,j+1,k)*An
              p=p+Phi(i,j,k-1)*Ab
              p=p+Phi(i,j,k+1)*At
              p=p-RHS(i,j,k)
              p=abs(-p +Ap*Phi(i,j,k))
              R=max(R,abs(p))
             end do
         end do
      end do

      call MPI_AllReduce(MPI_IN_PLACE, R, 1, MPI_RP, MPI_MAX, glob_comm, ie)

   end subroutine Res3D



    subroutine Res2D(nx,ny,Phi,RHS,&
                     Aw,Ae,As,An,&
                     R)

      implicit none
      intrinsic mod, abs, max

      integer,parameter:: KND=RP,DIRICHLET=1,NEUMANN=2,PERIODIC=3
      integer, parameter :: We=1,Ea=2,So=3,No=4

      integer,intent(in)   :: nx,ny
      real(KND),dimension(0:nx+1,0:ny+1),intent(inout)::Phi
      real(KND),dimension(1:nx,1:ny),intent(in)::RHS
      real(KND),intent(in) :: Aw,Ae
      real(KND),intent(in) :: As,An
      real(KND),intent(out) :: R
      integer i,j,k,l
      real(KND) :: p,Ap
      integer ie

      R=0
      Ap = Aw + Ae + As + An
      
      do j=1,ny
         do i=1,nx
            p=0
            p=p+Phi(i-1,j)*Aw
            p=p+Phi(i+1,j)*Ae
            p=p+Phi(i,j-1)*As
            p=p+Phi(i,j+1)*An
            p=p-RHS(i,j)
            p=abs(-p +Ap*Phi(i,j))
            R=max(R,abs(p))
           end do
       end do

      call MPI_AllReduce(MPI_IN_PLACE, R, 1, MPI_RP, MPI_MAX, glob_comm, ie)

   end subroutine Res2D





end module subs


















program testpoisson_MPI
  use iso_fortran_env
  use iso_c_binding
  use PoisFFT, PoisFFT_Solver1D => PoisFFT_Solver1D_DP, &
               PoisFFT_Solver2D => PoisFFT_Solver2D_DP, &
               PoisFFT_Solver3D => PoisFFT_Solver3D_DP
  use my_mpi
  use subs
  
  implicit none

  integer(c_intptr_t) :: ng(3) = [131, 123, 127]![21,32,25]!
  real(RP), dimension(:,:,:), allocatable :: Phi, RHS
  real(RP) :: dx,dy,dz, Ls(3)
  integer i,j,k
  real(RP) :: R,x,y,z
  integer(int64) :: t1,t2,trate
  type(PoisFFT_Solver3D) :: Solver3D
  type(PoisFFT_Solver2D) :: Solver2D
  type(PoisFFT_Solver1D) :: Solver1D
  character(len=50) :: ch50
  integer :: nx, ny, nz
  integer(c_intptr_t),dimension(3) :: nxyz,off,nxyz2,nsxyz2
  real(DRP) ::  S
  integer :: ie = 0
  character(200) :: fname

  integer :: required, provided

  required = MPI_THREAD_SERIALIZED

!$ if (.false.) then
    call MPI_Init_thread(required, provided, ie)
    provided = required
!$ else
!$  call MPI_Init_thread(required, provided, ie)
!$ end if

  if (ie/=0) call error_stop("Error initializing MPI.")

  glob_comm = MPI_COMM_WORLD

  myim = this_image()
  myrank = myim - 1

  if (myrank==0) master = .true.

  if (provided<required) then
    if (master) write(*,*) "------------------------------"
    if (master) write(*,*) "Error, the provided MPI threading support smaller than required!"
    if (master) write(*,*) "required:", required
    if (master) write(*,*) "provided:", provided
    if (master) write(*,*) "Trying to continue anyway, but a crash is likely and the results will be questionable."
    if (master) write(*,*) "------------------------------"
  end if
  
  if (RP == kind(1.)) then
    MPI_RP = MPI_REAL
  else if (RP == kind(1.D0)) then
    MPI_RP = MPI_DOUBLE_PRECISION
  end if
  
  call system_clock(count_rate=trate)
 
  nims = num_images()

  if (command_argument_count()>=1) then
    call get_command_argument(1,value=ch50)
    read(ch50,*,iostat=ie) npxyz
    if (ie/=0) then
      write(*,*) "The process grid should be provided as 'npx,npy,npz' where nxp,npy and npz are integers."
      stop
    end if
  else
    npxyz(1) = 1
    npxyz(2) = nint(sqrt(real(nims)))
    npxyz(3) = nims / npxyz(2)
    if (master)    write (*,*) "Trying to decompose in",npxyz,"process grid."
  end if

  if (command_argument_count()>=2) then
    call get_command_argument(2,value=ch50)
    read(ch50,*) ng
  end if

  nxims = npxyz(1)
  nyims = npxyz(2)
  nzims = npxyz(3)
  
  if (product(npxyz)/=nims) then
    if (master) then
      write(*,*) "Could not decompose the processes to N x N grid."
      write(*,*) "Try a perfect square for number of processes."
    end if
    call error_stop(25)
  end if

  call PoisFFT_InitMPIGrid(glob_comm, npxyz(3:2:-1), cart_comm, ie)
  if (ie/=0) call error_stop(30)
  
  call get_image_coords

  call PoisFFT_LocalGridSize(3,ng,cart_comm,nxyz,off,nxyz2,nsxyz2)
  if (any(nxyz/=nxyz2).or.any(off/=nsxyz2)) call error_stop(40)

  if (any(nxyz==0)) then
    write(*,*) "Process",pxyz,"has grid dimensions",nxyz,"."
    write(*,*) "Try different process grid distribution."
    call error_stop(45)
  end if
  
  nx = nxyz(1)
  ny = nxyz(2)
  nz = nxyz(3)
 
  Ls = [2*pi, 2*pi, 2*pi]
  dx = Ls(1)/ng(1)
  dy = Ls(2)/ng(2)
  dz = Ls(3)/ng(3)
  
  allocate(RHS(nx,ny,nz),stat=ie)
  if (ie/=0) call error_stop(50)

  allocate(Phi(0:nx+1,0:ny+1,0:nz+1))
  if (ie/=0) call error_stop(60)


  call random_seed(put=[(0,i=1,100)])

  do k=1,nz
   do j=1,ny
    do i=1,nx
     x=(i+off(1)-1._RP/2)*dx
     y=(j+off(2)-1._RP/2)*dy
     z=(k+off(3)-1._RP/2)*dz
     call random_number(RHS(i,j,k))
     call random_number(Phi(i,j,k))
!      RHS(i,j,k) = x * sin(y) / (abs(z)+0.1)
!      PHI(i,j,k) = x * sin(y) / (abs(z)+0.1)
    end do
   end do
  end do
  
  call MPI_AllReduce(sum(RHS), S, 1, MPI_RP, MPI_SUM, glob_comm, ie)
 
  RHS = RHS - S / product(ng)

  
  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "3D PNsNs"

  call compute3d([(PoisFFT_Periodic, i=1,2),(PoisFFT_NeumannStag, i=3,6)])


  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "3D PPNs"

  call compute3d([(PoisFFT_Periodic, i=1,4),(PoisFFT_NeumannStag, i=5,6)])


  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "3D staggered Dirichlet"

  call compute3D([(PoisFFT_DirichletStag, i=1,6)])


  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "3D staggered Neumann:"

  call compute3D([(PoisFFT_NeumannStag, i=1,6)])


  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "3D Periodic"

  call compute3D([(PoisFFT_Periodic, i=1,6)])
 
 
  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "2D Periodic"

  call compute2D([(PoisFFT_Periodic, i=1,4)])
  

  call MPI_Barrier(glob_comm,ie)
  
  if (master) write(*,*) "1D Periodic"

  call compute1D([(PoisFFT_Periodic, i=1,2)])


  call MPI_Barrier(glob_comm,ie)

  call save_vtk  
  
  call MPI_finalize(ie)
  
  
contains

  subroutine compute3D(BCs)
    integer, intent(in) :: BCs(6)
    
    Phi(1:nx,1:ny,1:nz) = RHS

    call exchange_boundaries_3D(glob_comm, Phi, nx, ny, nz, BCs)

    call Res3D(nx,ny,nz,Phi,RHS,&
                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
                 R)

    if (master) write (*,*) "R1:",R

    Solver3D = PoisFFT_Solver3D([nx,ny,nz],Ls,BCs,PoisFFT_FiniteDifference2, &
                                int(ng),int(off),cart_comm)

    do i=1,2
      call system_clock(count=t1)

      call Execute(Solver3D,Phi,RHS)

      call system_clock(count=t2)
      if (master) write(*,*) "solver cpu time", real(t2-t1)/real(trate)
    end do

    call Finalize(Solver3D)

    call exchange_boundaries_3D(glob_comm, Phi, nx, ny, nz, BCs)

    call Res3D(nx,ny,nz,Phi,RHS,&
                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),dz**(-2),dz**(-2),&
                 R)

    if (master) write (*,*) "R2:",R

    if (master) write(*,*) "--------"
  end subroutine

  
  subroutine compute2D(BCs)
    integer, intent(in) :: BCs(4)
    integer :: sub_comm, ie
    integer :: n
    
    call MPI_Cart_sub(cart_comm, [.false., .true.], sub_comm, ie)

    call MPI_AllReduce(sum(RHS(:,:,1)), S, 1, MPI_RP, MPI_SUM, sub_comm, ie)
 
    RHS(:,:,1) = RHS(:,:,1) - S / product(ng(1:2))

    call exchange_boundaries_2D(glob_comm, Phi(:,:,1), nx, ny, BCs)

    call Res2D(nx,ny,Phi(:,:,1),RHS(:,:,1),&
                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),&
                 R)

    if (master) write (*,*) "R1:",R

    Solver2D = PoisFFT_Solver2D([nx,ny],Ls(1:2),BCs,PoisFFT_FiniteDifference2, &
                                int(ng(1:2)),int(off(1:2)),sub_comm)

    do i=1,1
      call system_clock(count=t1)

      call Execute(Solver2D,Phi(:,:,i),RHS(:,:,i))

      call system_clock(count=t2)
      if (master) write(*,*) "solver cpu time", real(t2-t1)/real(trate)
    end do

    call Finalize(Solver2D)

    call exchange_boundaries_2D(glob_comm, Phi(:,:,1), nx, ny, BCs)

    call Res2D(nx,ny,Phi(:,:,1),RHS(:,:,1),&
                 dx**(-2),dx**(-2),dy**(-2),dy**(-2),&
                 R)

    if (master) write (*,*) "R2:",R

    if (master) write(*,*) "--------"
  end subroutine

  
  subroutine compute1D(BCs)
    integer, intent(in) :: BCs(2)
    integer :: sub_comm = -1, ie
    
    call MPI_Cart_sub(cart_comm, [.false., .true.], sub_comm, ie)

    call MPI_AllReduce(sum(RHS(1,:,1)), S, 1, MPI_RP, MPI_SUM, sub_comm, ie)
 
    RHS(1,:,1) = RHS(1,:,1) - S / ng(2)

!     call exchange_boundaries_1D(glob_comm, Phi(:,:,1), nx, BCs)

!     call Res1D(nx,ny,Phi(:,1,1),RHS(:,1,1),&
!                  dx**(-2),dx**(-2),&
!                  R)

    if (master) write (*,*) "R1:",R

    Solver1D = PoisFFT_Solver1D([ny],Ls(1:1),BCs,PoisFFT_FiniteDifference2, &
                                int(ng(2:2)),int(off(2:2)),sub_comm)

    do i=1,1
      call system_clock(count=t1)

      call Execute(Solver1D,Phi(1,:,i),RHS(1,:,i))

      call system_clock(count=t2)
      if (master) write(*,*) "solver cpu time", real(t2-t1)/real(trate)
    end do

    call Finalize(Solver1D)

!     call exchange_boundaries_1D(glob_comm, Phi(:,1,1), nx, ny, BCs)

!     call Res1D(nx,ny,Phi(:,1,1),RHS(:,1,1),&
!                  dx**(-2),dx**(-2),&
!                  R)

    if (master) write (*,*) "R2:",R

    if (master) write(*,*) "--------"
  end subroutine

  
  subroutine save_vtk
    use Endianness
    integer filetype, fh, unit
    integer(MPI_OFFSET_KIND) pos
    real(RP),allocatable :: buffer(:,:,:),buf(:)
    character :: lf=achar(10)
    character(70) :: str
    character(10) :: fm = '(*(1x,g0))'

    call GetEndianness
    
    if (master) then
      open(newunit=unit,file="out.vtk", &
           access='stream',status='replace',form="unformatted",action="write")

      write(unit) "# vtk DataFile Version 2.0",lf
      write(unit) "CLMM output file",lf
      write(unit) "BINARY",lf
      write(unit) "DATASET RECTILINEAR_GRID",lf
      str="DIMENSIONS"
      write(str(12:),fm) ng(1),ng(2),ng(3)
      write(unit) str,lf
      str="X_COORDINATES"
      write(str(15:),fm) ng(1),"double"
      write(unit) str,lf
      write(unit) BigEnd([((i-0.5)*dx,i=1,ng(1))]),lf
      str="Y_COORDINATES"
      write(str(15:),fm) ng(2),"double"
      write(unit) str,lf
      write(unit) BigEnd([((j-0.5)*dy,j=1,ng(2))]),lf
      str="Z_COORDINATES"
      write(str(15:),fm) ng(3),"double"
      write(unit) str,lf
      write(unit) BigEnd([((k-0.5)*dz,k=1,ng(3))]),lf
      str="POINT_DATA"
      write(str(12:),fm) product(ng)
      write(unit) str,lf
      write(unit) lf
      write(unit) "SCALARS Phi double",lf
      write(unit) "LOOKUP_TABLE default",lf
      close(unit)
    end if
    
    call MPI_Barrier(glob_comm,ie)
    
    call MPI_File_open(glob_comm,"out.vtk", MPI_MODE_APPEND + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ie)
    if (ie/=0) call error_stop("open")

    call MPI_Type_create_subarray(3, int(ng), int(nxyz), int(off), &
       MPI_ORDER_FORTRAN, MPI_RP, filetype, ie)
    if (ie/=0) call error_stop("create_subarray")
    
    call MPI_type_commit(filetype, ie)
    if (ie/=0) call error_stop("type_commit")
    
    call MPI_Barrier(glob_comm,ie)
    call MPI_File_get_position(fh, pos, ie)
    call MPI_Barrier(glob_comm,ie)
    
    call MPI_File_set_view(fh, pos, MPI_RP, filetype, "native", MPI_INFO_NULL, ie)
    if (ie/=0) call error_stop("set_view")
    
    allocate(buffer(1:nx,1:ny,1:nz))
    buffer = BigEnd(Phi(1:nx,1:ny,1:nz))
       
    call MPI_File_write_all(fh, buffer, nx*ny*nz, MPI_RP, MPI_STATUS_IGNORE, ie)
    if (ie/=0) call error_stop("write_all")
    
    call MPI_File_close(fh, ie)
    if (ie/=0) call error_stop("close")
    
  end subroutine

end program testpoisson_MPI
