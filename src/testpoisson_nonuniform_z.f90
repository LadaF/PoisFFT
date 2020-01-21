module Kinds
  use PoisFFT_Precisions
  use iso_c_binding
  use iso_fortran_env

  integer,parameter :: rp = DRP
  integer,parameter :: cp = DCP
  integer,parameter :: size_kind = c_size_t
end module

module Globals
  use Kinds
  
  implicit none

  real(rp), parameter :: pi = 3.141592653589793238462_rp!pi = 4*atan(1._rp)!
  real(rp), parameter :: Lx = 2*pi, Ly = Lx*1.1_rp, Lz = Lx/1.1_rp
  integer :: nx = 21, ny = 32, nz = 25
  real(rp) :: dx,dy,dz
  real(rp), dimension(:,:,:), allocatable :: Phi3D, RHS3D
  real(rp), dimension(:,:), allocatable :: Phi2D, RHS2D
  real(rp), dimension(:), allocatable :: Phi1D, RHS1D

  real(rp), dimension(:), allocatable :: z(:), z_u(:)
end module



module Residues
  use Kinds
  use Globals
  use PoisFFT
    
  implicit none
  
  integer, parameter :: We = 1,Ea = 2,So = 3,No = 4,Bo = 5,To = 6

contains


  
  
  subroutine BCond(Phi)
    real(rp), intent(inout) :: Phi(0:,0:,0:)
    Phi(0,:,:) = Phi(nx,:,:)
    Phi(nx+1,:,:) = Phi(1,:,:)
    Phi(:,0,:) = Phi(:,ny,:)
    Phi(:,ny+1,:) = Phi(:,1,:)
    Phi(:,:,0) = Phi(:,:,nz)
    Phi(:,:,nz+1) = Phi(:,:,1)
  end subroutine


  subroutine Res3D(Phi,RHS,&
                   Aw,Ae,As,An,z,z_u,&
                   Btype,R)
    !2nd order finite difference residuum

    real(rp), intent(inout) :: Phi(0:,0:,0:)
    real(rp), intent(in) :: RHS(:,:,:)
    real(rp),intent(in) :: Aw,Ae
    real(rp),intent(in) :: As,An
    real(rp),intent(in) :: z(0:),z_u(-1:)
    integer,intent(in) :: Btype(6)
    real(rp),intent(out) :: R
    integer i,j,k
    real(rp) :: p,Ap

    call BCond(Phi)
    
    R = 0

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          p = 0
          Ap = 0          
          if (i>1) then
                    p = p + Phi(i-1,j,k)*Aw
                    Ap = Ap + Aw
          elseif (Btype(We)==PoisFFT_Periodic) then
                    p = p + Phi(nx,j,k)*Aw
                    Ap = Ap + Aw
          elseif (Btype(We)==PoisFFT_DirichletStag) then
                    p = p - Phi(1,j,k)*Aw
                    Ap = Ap + Aw
          elseif (Btype(We)==PoisFFT_Dirichlet) then
                    Ap = Ap + Aw
          elseif (Btype(We)==PoisFFT_Neumann) then
                    p = p + Phi(2,j,k)*Aw
                    Ap = Ap + Aw
          end if
          if (i<nx) then
                    p = p + Phi(i+1,j,k)*Ae
                    Ap = Ap + Ae
          elseif (Btype(Ea)==PoisFFT_Periodic) then
                    p = p + Phi(1,j,k)*Ae
                    Ap = Ap + Ae
          elseif (Btype(Ea)==PoisFFT_DirichletStag) then
                    p = p - Phi(nx,j,k)*Ae
                    Ap = Ap + Ae
          elseif (Btype(Ea)==PoisFFT_Dirichlet) then
                    Ap = Ap + Ae
          elseif (Btype(Ea)==PoisFFT_Neumann) then
                    p = p + Phi(nx-1,j,k)*Ae
                    Ap = Ap + Ae
          end if
          if (j>1) then
                    p = p + Phi(i,j-1,k)*As
                    Ap = Ap + As
          elseif (Btype(So)==PoisFFT_Periodic) then
                    p = p + Phi(i,ny,k)*As
                    Ap = Ap + As
           elseif (Btype(So)==PoisFFT_DirichletStag) then
                    p = p - Phi(i,1,k)*As
                    Ap = Ap + As
          elseif (Btype(So)==PoisFFT_Dirichlet) then
                    Ap = Ap + As
          elseif (Btype(So)==PoisFFT_Neumann) then
                    p = p + Phi(i,2,k)*As
                    Ap = Ap + As
          end if
          if (j<ny) then
                    p = p + Phi(i,j+1,k)*An
                    Ap = Ap + An
          elseif (Btype(No)==PoisFFT_Periodic) then
                    p = p + Phi(i,1,k)*An
                    Ap = Ap + An
          elseif (Btype(No)==PoisFFT_DirichletStag) then
                    p = p - Phi(i,ny,k)*An
                    Ap = Ap + An
          elseif (Btype(No)==PoisFFT_Dirichlet) then
                    Ap = Ap + An
          elseif (Btype(No)==PoisFFT_Neumann) then
                    p = p + Phi(i,ny-1,k)*An
                    Ap = Ap + An
          end if
          if (k>1) then
                    p = p + Phi(i,j,k-1)/(z(k) - z(k-1))/(z_u(k)-z_u(k-1))
                    Ap = Ap + 1._rp/(z(k) - z(k-1))/(z_u(k)-z_u(k-1))
!           elseif (Btype(Bo)==PoisFFT_Periodic) then
!                     p = p + Phi(i,j,nz)*Ab
!                     Ap = Ap + Ab
!           elseif (Btype(Bo)==PoisFFT_DirichletStag) then
!                     p = p - Phi(i,j,1)*Ab
!                     Ap = Ap + Ab
!           elseif (Btype(Bo)==PoisFFT_Dirichlet) then
!                     Ap = Ap + Ab
!           elseif (Btype(Bo)==PoisFFT_Neumann) then
!                     p = p + Phi(i,j,2)*Ab
!                     Ap = Ap + Ab
          end if
          if (k<nz) then
                    p = p + Phi(i,j,k+1)/(z(k+1) - z(k))/(z_u(k)-z_u(k-1))
                    Ap = Ap + 1._rp/(z(k+1) - z(k))/(z_u(k)-z_u(k-1))
!           elseif (Btype(To)==PoisFFT_Periodic) then
!                     p = p + Phi(i,j,1)*At
!                     Ap = Ap + At
!           elseif (Btype(To)==PoisFFT_DirichletStag) then
!                     p = p - Phi(i,j,nz)*At
!                     Ap = Ap + At
!           elseif (Btype(To)==PoisFFT_Dirichlet) then
!                     Ap = Ap + At
!           elseif (Btype(To)==PoisFFT_Neumann) then
!                     p = p + Phi(i,j,nz-1)*At
!                     Ap = Ap + At
          end if

          p = p - RHS(i,j,k)
          p = abs(-p +Ap*Phi(i,j,k))
          R = max(R,abs(p))
        end do
      end do
    end do

  end subroutine Res3D
end module










module Tests
  use Kinds
  use Globals
  use Residues
  use PoisFFT, Solver_nonuniform_z => PoisFFT_Solver3D_nonuniform_z_DP
  
  implicit none
  
  character(20), parameter :: char_BCs(0:4) = ["Periodic            ", &
                                             "Dirichlet regular   ", &
                                             "Neumann regular     ", &
                                             "Dirichlet staggered ", &
                                             "Neumann staggered   "]
 
contains

  
  
  
  subroutine Test3D(BCs)
    integer, intent(in) :: BCs(6)
    integer :: i

    write(*,*) "----"
    write(*,'(1x,"3D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))

    call TestFD2_3D(BCs)
    write(*,*)
  end subroutine

  
  subroutine TestFD2_3D(BCs)
    integer, intent(in) :: BCs(6)
    real(rp) :: R
    
    call RunFD2_3D(BCs)
    
    call Res3D(Phi3D, RHS3D, &
               dx**(-2), dx**(-2), dy**(-2), dy**(-2), z, z_u, &
               BCs, R)
    if (R < int(nx, int64) * int(ny, int64) * int(nz, int64) * 1000 * epsilon(1._rp)) then
      write(*,*) "FD2 OK", R
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "FD2 residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_3D(BCs)
    type(Solver_nonuniform_z) :: Solver
    integer, intent(in) :: BCs(6)
    integer :: i
  
    Solver = Solver_nonuniform_z([nx,ny,nz], [Lx,Ly], z(1:nz), z_u(0:nz), &
                              BCs, &
                              approximation=2)                   

    call Execute(Solver, Phi3D(1:nx,1:ny,1:nz), RHS3D)

    call Finalize(Solver)

  end subroutine
  

  
end module
















program testpoisson
  use Kinds
  use Globals
  use Residues
  use Tests
  use PoisFFT

  
  implicit none

  integer i,j,k,niters
  real(rp) :: xp, yp, zp
  character(len = 12) :: arg
  real(RP) :: avg, p
  
  
  call random_seed(size=j)
  call random_seed(put=[(i+1,i=1,j)])
  
  if (command_argument_count()>=3) then
    call get_command_argument(1,value = arg)
    read(arg,'(i12)') nx
    call get_command_argument(2,value = arg)
    read(arg,'(i12)') ny
    call get_command_argument(3,value = arg)
    read(arg,'(i12)') nz
  end if
  
  write(*,*) "nx",nx
  write(*,*) "ny",ny
  write(*,*) "nz",nz
 
  dx = Lx/nx;dy = Ly/ny;dz = Lz/nz
  
  allocate(z(0:nz+1))
  allocate(z_u(-1:nz+1))
  z_u(0:nz) = [(Lz/nz*i, i = 0, nz)]
  do i = 1, nz-1
    call random_number(p)
    p = (p - 0.5_rp)*0.5
     z_u(i) = z_u(i) + min(abs(z_u(i+1)-z_u(i)),abs(z_u(i)-z_u(i-1)))*p
  end do
  z_u(-1) = z_u(0) - (z_u(1)-z_u(0))
  z_u(nz+1) = z_u(nz) + (z_u(nz)-z_u(nz-1))
  do i = 0, nz+1
    z(i) = (z_u(i)+z_u(i-1))/2
  end do

  
  allocate(RHS3D(nx,ny,nz))
  allocate(Phi3D(0:nx+1,0:ny+1,0:nz+1))


  avg = 0
  do k = 1,nz
   do j = 1,ny
    do i = 1,nx
     xp = (i-1._rp/2)*dx
     yp = (j-1._rp/2)*dy
     zp = z(k)
     call random_number(RHS3D(i,j,k))
     avg = avg + RHS3D(i,j,k)*dx*dy*(z_u(k)-z_u(k-1))
     call random_number(Phi3D(i,j,k))
    end do
   end do
  end do

  avg = avg / (Lx*Ly*Lz)
  RHS3D = RHS3D - avg


  dx = Lx / nx
  dy = Ly / ny
  call Test3D([(PoisFFT_PERIODIC, i = 1,4),(PoisFFT_NeumannStag, i = 5,6)])

  call Test3D([(PoisFFT_NeumannStag, i = 1,6)])


end program testpoisson
