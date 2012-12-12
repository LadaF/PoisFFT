module vtkarray
!   use iso_fortran_env, only: real32
  implicit none

  integer,parameter :: real32=selected_real_kind(p=6,r=37)

  contains

    subroutine VtkArraySimple(fname,A)
      character(len=*)                 :: fname
      real(real32),dimension(1:,1:,1:) :: A
      integer                          :: nx,ny,nz
      integer                          :: i,j,k
      character(len=40)                :: str
      
      nx=Ubound(A,1)
      ny=Ubound(A,2)
      nz=Ubound(A,3)

      open(11,file=fname)
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
      write (11,"(A)") trim(str)
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') nx,"float"
      write (11,"(A)") str
      write (11,*) (i, i=1,nx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') ny,"float"
      write (11,"(A)") trim(str)
      write (11,*) (i, i=1,ny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') nz,"float"
      write (11,"(A)") trim(str)
      write (11,*) (i, i=1,nz)
      str="POINT_DATA"
      write (str(12:),*) nx*ny*nz
      write (11,"(A)") trim(str)


      write (11,"(A)") "SCALARS array float"
      write (11,"(A)") "LOOKUP_TABLE default"
      do k=1,nz
       do j=1,ny
        do i=1,nx
          write (11,*) A(i,j,k)
        enddo
        enddo
      enddo
      write (11,*)

    end subroutine VtkArraySimple
end module vtkarray
