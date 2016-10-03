module Kinds
  use PoisFFT_Precisions
  use iso_c_binding
  use iso_fortran_env

  integer,parameter :: rp = DRP
  integer,parameter :: size_kind = c_size_t
end module

module Globals
  use Kinds
  
  implicit none

  real(rp), parameter :: pi = 3.141592653589793238462_rp!pi = 4*atan(1._rp)!
  real(rp), parameter :: Lx = 2*pi, Ly = Lx*1.1, Lz = Lx/1.1
  integer :: nx = 21, ny = 32, nz = 25
  real(rp) :: dx,dy,dz
  real(rp), dimension(:,:,:), allocatable :: Phi3D, RHS3D
  real(rp), dimension(:,:), allocatable :: Phi2D, RHS2D
  real(rp), dimension(:), allocatable :: Phi1D, RHS1D

end module



module Residues
  use Kinds
  use Globals
  use PoisFFT
  
  implicit none
  

contains

  subroutine ResExact1D_Dir(Phi, R)
    real(rp),intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f
    xx(i) = dx*i
!       f(x) = (x**3 - x*(nx*dx)**2)/6
    f(x) = -Lx**2/(5*pi)**2 * sin(5*pi*x/Lx)
!     f(x) = -(x/(3*Lx**2))*(Lx**3-2*Lx*x**2+x**3)
    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
!         print *, x, f(x), Phi(i)
      R = R + p**2
    end do
    R = sqrt(R)/nx
  end subroutine

  subroutine ResExact1D_DirStag(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f
    xx(i) = dx/2 + dx*(i-1)
!       f(x) = (x**3 - x*(nx*dx)**2)/6
    f(x) = -Lx**2/(5*pi)**2 * sin(5*pi*x/Lx)
!     f(x) = -(x/(3*Lx**2))*(Lx**3-2*Lx*x**2+x**3)
    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
      R = R + p**2
    end do
    R = sqrt(R)/nx
  end subroutine

  subroutine ResExact1D_Neum(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f,g
!       xx(i) = dx/2 + dx*(i-1) - 0.5*Lx
!       g(x) = (x**3/6 - x*(nx*dx)**2/8)
!       g(x) = (5._rp/16._rp)*Lx**2*x+x**5/(5*Lx**2)-x**3/2
    xx(i) = dx*(i-1)
    g(x) = -(Lx/(6*pi))**2 * cos(6*pi*(x/Lx))
    f(x) = g(x) + (Phi(1)-g(xx(1)))

    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
!         print *, x, f(x), Phi(i)
      R = R + p**2
    end do
    R = sqrt(R)/nx
  end subroutine


  subroutine ResExact1D_NeumStag(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f,g
!       xx(i) = dx/2 + dx*(i-1) - 0.5*Lx
!       g(x) = (x**3/6 - x*(nx*dx)**2/8)
!       g(x) = (5._rp/16._rp)*Lx**2*x+x**5/(5*Lx**2)-x**3/2
    xx(i) = dx/2 + dx*(i-1)
    g(x) = -(Lx/(6*pi))**2 * cos(6*pi*(x/Lx))
    f(x) = g(x) + (Phi(1)-g(xx(1)))

    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
      R = R + p**2
    end do
    R = sqrt(R)/nx
  end subroutine


  subroutine ResExact1D_per(Phi, R)
    real(rp), intent(in) :: Phi(0:)
    real(rp), intent(out) :: R
    integer :: i
    real(rp) :: x, p
    real(rp) :: xx,f,g
    xx(i) = dx*(i-1)
    g(x) = -(Lx/(6*pi))**2 * cos(6*pi*(x/Lx))
    f(x) = g(x) + (Phi(1)-g(xx(1)))

    R = 0
    do i = 1, nx
      x = xx(i)
      p = abs( Phi(i) - f(x) )
!         print *, x, f(x), Phi(i)
      R = R + p**2
    end do
    R = sqrt(R)/nx
  end subroutine
  
  
  
  
  
  
  
  
  
  subroutine ResExact3D_Dir(Phi, R)
    real(rp),intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx*i
    yy(j) = dy*j
    zz(k) = dz*k

    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (7*pi)**2/Lz**2 ) * &
                  sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)

    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R)/nx
  end subroutine

  subroutine ResExact3D_DirStag(Phi, R)
    real(rp), intent(in) :: Phi(0:,0:,0:)
    real(rp), intent(out) :: R
    integer :: i, j, k
    real(rp) :: x, y, z, p
    real(rp) :: xx, yy, zz, f
    
    xx(i) = dx/2 + dx*(i-1)
    yy(j) = dy/2 + dy*(j-1)
    zz(k) = dz/2 + dz*(k-1)
    
    f(x,y,z) = - 1 / ( (3*pi)**2/Lx**2 + (5*pi)**2/Ly**2 + (7*pi)**2/Lz**2 ) * &
                  sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)
                  
    R = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x = xx(i)
          y = yy(j)
          z = zz(k)
          p = abs( Phi(i, j, k) - f(x,y,z) )
          R = R + p**2
        end do
      end do
    end do
    R = sqrt(R)/nx
  end subroutine


  
  
  
  
  
  
  
  
  
  
  
  


  subroutine Res1D(Phi,RHS,&
                   Aw,Ae,&
                   Btype,R)
    !2nd order finite difference residuum

    integer, parameter :: Ea = 1,We = 2

    real(rp), intent(inout) :: Phi(0:)
    real(rp), intent(in) :: RHS(:)
    real(rp),intent(in) :: Aw,Ae
    integer,intent(in) :: Btype(2)
    real(rp),intent(out) :: R
    integer i
    real(rp) :: p,Ap

    R = 0

    do i = 1, nx
      p = 0
      Ap = 0
      if (i>1) then
                p = p + Phi(i-1)*Aw
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_Periodic) then
                p = p + Phi(nx)*Aw
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_DirichletStag) then
                p = p - Phi(1)*Aw
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_Dirichlet) then
                Ap = Ap + Aw
      elseif (Btype(We)==PoisFFT_Neumann) then
                p = p + Phi(2)*Aw
                Ap = Ap + Aw
      end if
      if (i<nx) then
                p = p + Phi(i+1)*Ae
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_Periodic) then
                p = p + Phi(1)*Ae
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_DirichletStag) then
                p = p - Phi(nx)*Ae
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_Dirichlet) then
                Ap = Ap + Ae
      elseif (Btype(Ea)==PoisFFT_Neumann) then
                p = p + Phi(nx-1)*Ae
                Ap = Ap + Ae
      end if

      p = p - RHS(i)
      p = abs(-p +Ap*Phi(i))
      R = max(R,abs(p))
    end do

  end subroutine Res1D










  subroutine Res2D(Phi,RHS,&
                   Aw,Ae,As,An,&
                   Btype,R)
    !2nd order finite difference residuum

    integer, parameter :: Ea = 1,We = 2,So = 3,No = 4

    real(rp), intent(inout) :: Phi(0:,0:)
    real(rp), intent(in) :: RHS(:,:)
    real(rp),intent(in) :: Aw,Ae
    real(rp),intent(in) :: As,An
    integer,intent(in) :: Btype(4)
    real(rp),intent(out) :: R
    integer i,j
    real(rp) :: p,Ap

    R = 0

    do j = 1, ny
      do i = 1, nx
        p = 0
        Ap = 0
        if (i>1) then
                  p = p + Phi(i-1,j)*Aw
                  Ap = Ap + Aw
        elseif (Btype(We)==PoisFFT_Periodic) then
                  p = p + Phi(nx,j)*Aw
                  Ap = Ap + Aw
        elseif (Btype(We)==PoisFFT_DirichletStag) then
                  p = p - Phi(1,j)*Aw
                  Ap = Ap + Aw
        elseif (Btype(We)==PoisFFT_Dirichlet) then
                  Ap = Ap + Aw
        elseif (Btype(We)==PoisFFT_Neumann) then
                  p = p + Phi(2,j)*Aw
                  Ap = Ap + Aw
        end if
        if (i<nx) then
                  p = p + Phi(i+1,j)*Ae
                  Ap = Ap + Ae
        elseif (Btype(Ea)==PoisFFT_Periodic) then
                  p = p + Phi(1,j)*Ae
                  Ap = Ap + Ae
        elseif (Btype(Ea)==PoisFFT_DirichletStag) then
                  p = p - Phi(nx,j)*Ae
                  Ap = Ap + Ae
        elseif (Btype(Ea)==PoisFFT_Dirichlet) then
                  Ap = Ap + Ae
        elseif (Btype(Ea)==PoisFFT_Neumann) then
                  p = p + Phi(nx-1,j)*Ae
                  Ap = Ap + Ae
        end if
        if (j>1) then
                  p = p + Phi(i,j-1)*As
                  Ap = Ap + As
        elseif (Btype(So)==PoisFFT_Periodic) then
                  p = p + Phi(i,ny)*As
                  Ap = Ap + As
        elseif (Btype(So)==PoisFFT_DirichletStag) then
                  p = p - Phi(i,1)*As
                  Ap = Ap + As
        elseif (Btype(So)==PoisFFT_Dirichlet) then
                  Ap = Ap + As
        elseif (Btype(So)==PoisFFT_Neumann) then
                  p = p + Phi(i,2)*As
                  Ap = Ap + As
        end if
        if (j<ny) then
                  p = p + Phi(i,j+1)*An
                  Ap = Ap + An
        elseif (Btype(No)==PoisFFT_Periodic) then
                  p = p + Phi(i,1)*An
                  Ap = Ap + An
        elseif (Btype(No)==PoisFFT_DirichletStag) then
                  p = p - Phi(i,ny)*An
                  Ap = Ap + An
        elseif (Btype(So)==PoisFFT_Dirichlet) then
                  Ap = Ap + An
        elseif (Btype(So)==PoisFFT_Neumann) then
                  p = p + Phi(i,ny-1)*An
                  Ap = Ap + An
        end if

        p = p - RHS(i,j)
        p = abs(-p +Ap*Phi(i,j))
        R = max(R,abs(p))
      end do
    end do

  end subroutine Res2D









  subroutine Res3D(Phi,RHS,&
                   Aw,Ae,As,An,Ab,At,&
                   Btype,R)
    !2nd order finite difference residuum

    integer, parameter :: Ea = 1,We = 2,So = 3,No = 4,Bo = 5,To = 6

    real(rp), intent(inout) :: Phi(0:,0:,0:)
    real(rp), intent(in) :: RHS(:,:,:)
    real(rp),intent(in) :: Aw,Ae
    real(rp),intent(in) :: As,An
    real(rp),intent(in) :: Ab,At
    integer,intent(in) :: Btype(6)
    real(rp),intent(out) :: R
    integer i,j,k
    real(rp) :: p,Ap

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
                    p = p + Phi(i,j,k-1)*Ab
                    Ap = Ap + Ab
          elseif (Btype(Bo)==PoisFFT_Periodic) then
                    p = p + Phi(i,j,nz)*Ab
                    Ap = Ap + Ab
           elseif (Btype(Bo)==PoisFFT_DirichletStag) then
                    p = p - Phi(i,j,1)*Ab
                    Ap = Ap + Ab
          elseif (Btype(Bo)==PoisFFT_Dirichlet) then
                    Ap = Ap + Ab
          elseif (Btype(Bo)==PoisFFT_Neumann) then
                    p = p + Phi(i,j,2)*Ab
                    Ap = Ap + Ab
          end if
          if (k<nz) then
                    p = p + Phi(i,j,k+1)*At
                    Ap = Ap + At
          elseif (Btype(To)==PoisFFT_Periodic) then
                    p = p + Phi(i,j,1)*At
                    Ap = Ap + At
          elseif (Btype(To)==PoisFFT_DirichletStag) then
                    p = p - Phi(i,j,nz)*At
                    Ap = Ap + At
          elseif (Btype(To)==PoisFFT_Dirichlet) then
                    Ap = Ap + At
          elseif (Btype(To)==PoisFFT_Neumann) then
                    p = p + Phi(i,j,nz-1)*At
                    Ap = Ap + At
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
  use PoisFFT, PoisFFT_Solver1D => PoisFFT_Solver1D_DP, &
               PoisFFT_Solver2D => PoisFFT_Solver2D_DP, &
               PoisFFT_Solver3D => PoisFFT_Solver3D_DP  
  
  implicit none
  
  character(20), parameter :: char_BCs(0:4) = ["Periodic            ", &
                                             "Dirichlet regular   ", &
                                             "Neumann regular     ", &
                                             "Dirichlet staggered ", &
                                             "Neumann staggered   "]
 
contains
  
  subroutine Test1D(BCs, test_proc)
    integer, intent(in) :: BCs(2)
    integer :: i
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    write(*,*) "----"
    write(*,'(1x,"1D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))
    
    call TestSpectral1D(BCs, test_proc)
    call TestFD2_1D(BCs)
    write(*,*)
  end subroutine
  
  subroutine TestSpectral1D(BCs, test_proc)
    integer, intent(in) :: BCs(2)
    real(rp) :: R
    
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    call RunSpectral1D(BCs)
    
    call test_proc(Phi1D, R)
    
    if (R < nx * 10 * epsilon(1._rp)) then
      write(*,*) "Spectral OK"
    else
      write(*,*) "Spectral FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunSpectral1D(BCs)
    type(PoisFFT_Solver1D) :: Solver
    integer, intent(in) :: BCs(2)
    
    Solver = PoisFFT_Solver1D([nx],[Lx],BCs)

    call Execute(Solver, Phi1D, RHS1D)

    call Finalize(Solver)
  end subroutine

  
  subroutine TestFD2_1D(BCs)
    integer, intent(in) :: BCs(2)
    real(rp) :: R
    
    call RunFD2_1D(BCs)
    
    call Res1D(Phi1D, RHS1D, &
               dx**(-2), dx**(-2), &
               BCs, R)
    
    if (R < nx * 10 * epsilon(1._rp)) then
      write(*,*) "FD2 OK"
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "Finite difference 2 residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_1D(BCs)
    type(PoisFFT_Solver1D) :: Solver
    integer, intent(in) :: BCs(2)
    
    Solver = PoisFFT_Solver1D([nx],[Lx],BCs, approximation=2)

    call Execute(Solver, Phi1D, RHS1D)

    call Finalize(Solver)
  end subroutine
  
  
  
  
  
  
  
  
  
  subroutine Test2D(BCs)!, test_pro
    integer, intent(in) :: BCs(4)
    integer :: i
!     interface
!       subroutine  test_proc(Phi, R)
!         import
!         real(rp),intent(in) :: Phi(:,:)
!         real(rp), intent(out) :: R
!       end subroutine
!     end interface
    
    write(*,*) "----"
    write(*,'(1x,"2D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))
    
!     call TestSpectral2D(BCs, test_proc)
    call TestFD2_2D(BCs)
    write(*,*)
  end subroutine
  
  subroutine TestSpectral2D(BCs, test_proc)
    integer, intent(in) :: BCs(4)
    real(rp) :: R
    
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:,:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    call RunSpectral2D(BCs)
    
    call test_proc(Phi2D, R)
    
    if (R < int(nx, int64) * int(ny, int64) * 10 * epsilon(1._rp)) then
      write(*,*) "Spectral OK"
    else
      write(*,*) "Spectral FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunSpectral2D(BCs)
    type(PoisFFT_Solver2D) :: Solver
    integer, intent(in) :: BCs(4)
    
    Solver = PoisFFT_Solver2D([nx,ny],[Lx,Ly],BCs)

    call Execute(Solver, Phi2D, RHS2D)

    call Finalize(Solver)
  end subroutine

  
  subroutine TestFD2_2D(BCs)
    integer, intent(in) :: BCs(4)
    real(rp) :: R
    
    call RunFD2_2D(BCs)
    
    call Res2D(Phi2D, RHS2D, &
               dx**(-2), dx**(-2), dy**(-2), dy**(-2), &
               BCs, R)
    if (R < int(nx, int64) * int(ny, int64) * 10 * epsilon(1._rp)) then
      write(*,*) "FD2 OK"
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_2D(BCs)
    type(PoisFFT_Solver2D) :: Solver
    integer, intent(in) :: BCs(4)
    
    Solver = PoisFFT_Solver2D([nx,ny],[Lx,Ly],BCs, approximation=2)

    call Execute(Solver, Phi2D, RHS2D)

    call Finalize(Solver)
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine Test3D(BCs, test_proc)
    integer, intent(in) :: BCs(6)
    integer :: i
    procedure(ResExact3D_Dir),optional :: test_proc
    
    write(*,*) "----"
    write(*,'(1x,"3D  |",*(a,"|"))') (trim(char_BCs(BCs(i))), i=1,size(BCs))
    
    if (present(test_proc)) call TestSpectral3D(BCs, test_proc)
    call TestFD2_3D(BCs)
    write(*,*)
  end subroutine
  
  subroutine TestSpectral3D(BCs, test_proc)
    integer, intent(in) :: BCs(6)
    real(rp) :: R
    
    interface
      subroutine  test_proc(Phi, R)
        import
        real(rp),intent(in) :: Phi(:,:,:)
        real(rp), intent(out) :: R
      end subroutine
    end interface
    
    call RunSpectral3D(BCs)
    
    call test_proc(Phi3D, R)
    
    if (R < int(nx, int64) * int(ny, int64) * int(nz, int64) * 10 * epsilon(1._rp)) then
      write(*,*) "Spectral OK"
    else
      write(*,*) "Spectral FAIL"
      write(*,*) "Spectral residuum:",R
    end if
  end subroutine
  
  subroutine RunSpectral3D(BCs)
    type(PoisFFT_Solver3D) :: Solver
    integer, intent(in) :: BCs(6)
    
    Solver = PoisFFT_Solver3D([nx,ny,nz],[Lx,Ly,Lz],BCs)

    call Execute(Solver, Phi3D, RHS3D)

    call Finalize(Solver)
  end subroutine

  
  subroutine TestFD2_3D(BCs)
    integer, intent(in) :: BCs(6)
    real(rp) :: R
    
    call RunFD2_3D(BCs)
    
    call Res3D(Phi3D, RHS3D, &
               dx**(-2), dx**(-2), dy**(-2), dy**(-2), dz**(-2), dz**(-2), &
               BCs, R)
    if (R < int(nx, int64) * int(ny, int64) * int(nz, int64) * 10 * epsilon(1._rp)) then
      write(*,*) "FD2 OK"
    else
      write(*,*) "FD2 FAIL"
      write(*,*) "Finite difference 2 residuum:",R
    end if
  end subroutine
  
  subroutine RunFD2_3D(BCs)
    type(PoisFFT_Solver3D) :: Solver
    integer, intent(in) :: BCs(6)
    
    Solver = PoisFFT_Solver3D([nx,ny,nz],[Lx,Ly,Lz],BCs, approximation=2)

    call Execute(Solver, Phi3D, RHS3D)

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
  real(rp) :: x,y,z
  character(len = 12) :: arg
  real(RP) :: avg, p

  
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
 
  dx = Lx/nx;dy = Lx/ny;dz = Lx/nz

  allocate(RHS3D(nx,ny,nz))
  allocate(Phi3D(0:nx+1,0:ny+1,0:nz+1))

  allocate(RHS2D(nx,ny))
  allocate(Phi2D(0:nx+1,0:ny+1))

  allocate(RHS1D(nx))
  allocate(Phi1D(0:nx+1))



  do k = 1,nz
   do j = 1,ny
    do i = 1,nx
     x=(i-1._rp/2)*dx
     y=(j-1._rp/2)*dy
     z=(k-1._rp/2)*dz
     call random_number(RHS3D(i,j,k))
     call random_number(Phi3D(i,j,k))
    end do
   end do
  end do
  RHS3D = RHS3D-sum(RHS3D)/size(RHS3D)




  do j = 1,ny
    do i = 1,nx
     x=(i-1._rp/2)*dx
     y=(j-1._rp/2)*dy
     z=(k-1._rp/2)*dz
     call random_number(RHS2D(i,j))
     call random_number(Phi2D(i,j))
    end do
  end do
  RHS2D = RHS2D-sum(RHS2D)/size(RHS2D)



  do i = 1,nx
    x=(i-1._rp/2)*dx
    call random_number(RHS1D(i))
    call random_number(Phi1D(i))
  end do
  RHS1D = RHS1D-sum(RHS1D)/size(RHS1D)


  dx = Lx / nx
  do i = 1,nx
    x = dx*(i-1)
    RHS1D(i) = cos(3*2*pi*(x/Lx)) !+ cos(13*2*pi*(x/Lx))
  end do
  call Test1D([(PoisFFT_PERIODIC, i = 1,2)], ResExact1D_Per)

  
  dx = Lx / (nx+1)
!   RHS1D = [(dx/2 + dx*(i-1), i = 1,nx)]
  RHS1D = [(sin(5*pi*(dx*(i))/Lx), i = 1,nx)]
!   RHS1D = [(  -(dx*i-Lx/2)**2/(Lx/2)**2+1, i = 1,nx)]
  call Test1D([(PoisFFT_Dirichlet, i = 1,2)], ResExact1D_Dir)


  dx = Lx / nx
  RHS1D = [(sin(5*pi*(dx*(i-0.5_rp))/Lx), i = 1,nx)]
!   RHS1D = [(  -(dx/2 + dx*(i-1)-Lx/2)**2/(Lx/2)**2+1, i = 1,nx)]
  call Test1D([(PoisFFT_DirichletStag, i = 1,2)], ResExact1D_DirStag)


  dx = Lx / (nx-1)
  do i = 1,nx
    x = dx*(i-1)
    RHS1D(i) = cos(3*2*pi*(x/Lx)) !+ cos(13*2*pi*(x/Lx))
  end do
  RHS1D = RHS1D - (sum(RHS1D(2:nx-1))+RHS1D(1)/2+RHS1D(nx)/2)/(size(RHS1D,kind=size_kind)-1)
  call Test1D([(PoisFFT_Neumann, i = 1,2)], ResExact1D_Neum)


  dx = Lx / nx
  do i = 1,nx
    x = dx/2 + dx*(i-1)
    RHS1D(i) = cos(3*2*pi*(x/Lx)) !+ cos(13*2*pi*(x/Lx))
  end do
  RHS1D = RHS1D - sum(RHS1D)/size(RHS1D,kind=size_kind)
  call Test1D([(PoisFFT_NeumannStag, i = 1,2)], ResExact1D_NeumStag)

  



  dx = Lx / nx
  dy = Ly / ny
  call Test2D([(PoisFFT_Periodic, i = 1,4)])

  dx = Lx / (nx+1)
  dy = Ly / (ny+1)
  call Test2D([(PoisFFT_Dirichlet, i = 1,4)])

  dx = Lx / nx
  dy = Ly / ny
  call Test2D([(PoisFFT_DirichletStag, i = 1,4)])

  dx = Lx / (nx-1)
  dy = Ly / (ny-1)
  avg = 0
  do j = 1,ny
    do i = 1,nx
      p = RHS2D(i,j)
      if (i==1.or.i==nx) p = p / 2
      if (j==1.or.j==ny) p = p / 2
      avg = avg + p
    end do
  end do
  RHS2D = RHS2D - avg/(real(nx-1,rp)*real(ny-1,rp))
  call Test2D([(PoisFFT_Neumann, i = 1,4)])

  dx = Lx / nx
  dy = Ly / ny
  RHS2D = RHS2D - sum(RHS2D)/(size(RHS2D,kind=size_kind))
  call Test2D([(PoisFFT_NeumannStag, i = 1,4)])








  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  call Test3D([(PoisFFT_Periodic, i = 1,6)])


  dx = Lx / (nx+1)
  dy = Ly / (ny+1)
  dz = Lz / (nz+1)
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*i; y = dy*j; z = dz*k
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_Dirichlet, i = 1,6)], ResExact3D_Dir)


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        x = dx*(i-0.5_rp); y = dy*(j-0.5_rp); z = dz*(k-0.5_rp)
        RHS3D(i,j,k) = sin(3*pi*x/Lx) * sin(5*pi*y/Ly) *sin(7*pi*z/Lz)
      end do
    end do
  end do
  call Test3D([(PoisFFT_DirichletStag, i = 1,6)], ResExact3D_DirStag)


  dx = Lx / (nx-1)
  dy = Ly / (ny-1)
  dz = Lz / (nz-1)
  avg = 0
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        p = RHS3D(i,j,k)
        if (i==1.or.i==nx) p = p / 2
        if (j==1.or.j==ny) p = p / 2
        if (k==1.or.k==nz) p = p / 2
        avg = avg + p
      end do
    end do
  end do
  RHS3D = RHS3D - avg/(real(nx-1,rp)*real(ny-1,rp)*real(nz-1,rp))
  call Test3D([(PoisFFT_Neumann, i = 1,6)])


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  RHS3D = RHS3D - sum(RHS3D)/(size(RHS3D,kind=size_kind))
  call Test3D([(PoisFFT_NeumannStag, i = 1,6)])


  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  call Test3D([(PoisFFT_PERIODIC, i = 1,4),(PoisFFT_NeumannStag, i = 5,6)])

  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  call Test3D([(PoisFFT_PERIODIC, i = 1,2),(PoisFFT_NeumannStag, i = 3,4),(PoisFFT_PERIODIC, i = 5,6)])

  dx = Lx / nx
  dy = Ly / ny
  dz = Lz / nz
  call Test3D([(PoisFFT_PERIODIC, i = 1,2),(PoisFFT_NeumannStag, i = 3,6)])

end program testpoisson
