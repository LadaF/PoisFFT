module Endianness

  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd, SwapB, littleendian

  logical, save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd16
    module procedure BigEnd32
    module procedure BigEnd32_a1
    module procedure BigEnd32_a2
    module procedure BigEnd32_a3
    module procedure BigEnd32_a4
    module procedure BigEnd64
    module procedure BigEnd64_a1
    module procedure BigEnd64_a2
    module procedure BigEnd64_a3
    module procedure BigEnd64_a4
  end interface

  interface SwapB
    module procedure SwapB32
    module procedure SwapB64
  end interface

  contains

    subroutine GetEndianness
      character(4) :: bytes !may not work on some processors

      bytes = transfer(1_int32,bytes)
      if (ichar(bytes(4:4))==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

   
    elemental function BigEnd16(x) result(res)
      integer(int16) :: res
      integer(int16),intent(in)::x
      character(2) :: bytes
      
      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes)
        res = ichar(bytes(2:2),int16)
        res = ior( ishft(ichar(bytes(1:1),int16),8), res )
      endif
    end function
    
    function BigEnd32(x) result(res)
      real(real32),intent(in) :: x
      real(real32) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a1(x) result(res)
      real(real32),intent(in) :: x(:)
      real(real32) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a2(x) result(res)
      real(real32),intent(in) :: x(:,:)
      real(real32) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a3(x) result(res)
      real(real32),intent(in) :: x(:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a4(x) result(res)
      real(real32),intent(in) :: x(:,:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64(x) result(res)
      real(real64),intent(in) :: x
      real(real64) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a1(x) result(res)
      real(real64),intent(in) :: x(:)
      real(real64) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a2(x) result(res)
      real(real64),intent(in) :: x(:,:)
      real(real64) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a3(x) result(res)
      real(real64),intent(in) :: x(:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a4(x) result(res)
      real(real64),intent(in) :: x(:,:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    elemental function SwapB32(x) result(res)
      real(real32) :: res
      real(real32),intent(in) :: x
      character(4) :: bytes
      integer(int32) :: t
      real(real32) :: rbytes, rt
      !nicer looking TRANSFER is problematic to optimize by ifort
      ! gfortran would be OK with that
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(4:4),int32)

      t = ior( ishftc(ichar(bytes(3:3),int32),8),  t )

      t = ior( ishftc(ichar(bytes(2:2),int32),16), t )

      t = ior( ishftc(ichar(bytes(1:1),int32),24), t )

      res = rt
        
    end function
    
    elemental function SwapB64(x) result(res)
      real(real64) :: res
      real(real64),intent(in) :: x
      character(8) :: bytes
      integer(int64) :: t
      real(real64) :: rbytes, rt
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(8:8),int64)

      t = ior( ishftc(ichar(bytes(7:7),int64),8),  t )

      t = ior( ishftc(ichar(bytes(6:6),int64),16), t )

      t = ior( ishftc(ichar(bytes(5:5),int64),24), t )

      t = ior( ishftc(ichar(bytes(4:4),int64),32), t )

      t = ior( ishftc(ichar(bytes(3:3),int64),40), t )

      t = ior( ishftc(ichar(bytes(2:2),int64),48), t )

      t = ior( ishftc(ichar(bytes(1:1),int64),56), t )

      res = rt

    end function

end module Endianness

module vtkarray
 !Simple module to output arrays for visualization. No physical coordinates are used, only the position in the array.
 !Mostly only for debugging.

  use iso_fortran_env, only: real32, real64
  use Endianness, only: BigEnd

  implicit none

  interface VtkArrayAscii
    module procedure SVtkArrayAscii
    module procedure DVtkArrayAscii
  end interface

  interface VtkArrayBin
    module procedure SVtkArrayBin
    module procedure DVtkArrayBin
  end interface
  
  character, parameter :: lf = achar(10)

contains

  subroutine SVtkArrayAscii(fname, A)
    character(len=*), intent(in)  :: fname
    real(real32), intent(in)      :: A(:,:,:)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    character(len=40)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)

    open(newunit=unit,file=fname,status="replace",action="write")
    write(unit,"(A)") "# vtk DataFile Version 2.0"
    write(unit,"(A)") "CLMM output file"
    write(unit,"(A)") "ASCII"
    write(unit,"(A)") "DATASET RECTILINEAR_GRID"
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit,"(A)") trim(str)
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit,"(A)") str
    write(unit,*) (i, i=1,nx)
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,ny)
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,nz)
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit,"(A)") trim(str)


    write(unit,"(A)") "SCALARS array float"
    write(unit,"(A)") "LOOKUP_TABLE default"

    write(unit,'(*(g0,/))') A(1:nx,1:ny,1:nz)

    write(unit,*)
    close(unit)

  end subroutine SVtkArrayAscii

  subroutine DVtkArrayAscii(fname, A)
    character(len=*), intent(in)  :: fname
    real(real64), intent(in)      :: A(:,:,:)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    character(len=40)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)

    open(newunit=unit,file=fname,status="replace",action="write")
    write(unit,"(A)") "# vtk DataFile Version 2.0"
    write(unit,"(A)") "CLMM output file"
    write(unit,"(A)") "ASCII"
    write(unit,"(A)") "DATASET RECTILINEAR_GRID"
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit,"(A)") trim(str)
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"double"
    write(unit,"(A)") str
    write(unit,*) (i, i=1,nx)
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"double"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,ny)
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"double"
    write(unit,"(A)") trim(str)
    write(unit,*) (i, i=1,nz)
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit,"(A)") trim(str)


    write(unit,"(A)") "SCALARS array double"
    write(unit,"(A)") "LOOKUP_TABLE default"

    write(unit,'(*(g0,/))') A(1:nx,1:ny,1:nz)

    write(unit,*)
    close(unit)

  end subroutine DVtkArrayAscii
  
  subroutine SVtkArrayBin(fname, A, offsets)
    character(len=*), intent(in)  :: fname
    real(real32), intent(in)      :: A(:,:,:)
    integer, optional, intent(in) :: offsets(3)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    integer                       :: offs(3)
    character(len=70)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)
    
    if (present(offsets)) then
      offs = offsets
    else
      offs = 0
    end if

    open(newunit=unit,file=fname,access='stream',form='unformatted',status="replace",action="write")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(1), i=1,nx)], real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(2), i=1,ny)], real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(3), i=1,nz)], real32)), lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str, lf


    write(unit) "SCALARS array float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) Bigend(A(1:nx,1:ny,1:nz)), lf

    close(unit)

  end subroutine SVtkArrayBin

  subroutine DVtkArrayBin(fname, A, offsets)
    character(len=*), intent(in)  :: fname
    real(real64), intent(in)      :: A(:,:,:)
    integer, optional, intent(in) :: offsets(3)
    integer                       :: nx,ny,nz
    integer                       :: i
    integer                       :: unit
    integer                       :: offs(3)
    character(len=70)             :: str

    nx = ubound(A,1)
    ny = ubound(A,2)
    nz = ubound(A,3)
    
    if (present(offsets)) then
      offs = offsets
    else
      offs = 0
    end if

    open(newunit=unit,file=fname,access='stream',form='unformatted',status="replace",action="write")
    write(unit) "# vtk DataFile Version 2.0", lf
    write(unit) "CLMM output file", lf
    write(unit) "BINARY", lf
    write(unit) "DATASET RECTILINEAR_GRID", lf
    str="DIMENSIONS"
    write(str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
    write(unit) str, lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(1), i=1,nx)], real32)), lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(2), i=1,ny)], real32)), lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str, lf
    write(unit) BigEnd(real([(i+offs(3), i=1,nz)], real32)), lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str, lf


    write(unit) "SCALARS array float", lf
    write(unit) "LOOKUP_TABLE default", lf

    write(unit) BigEnd(real(A(1:nx,1:ny,1:nz), real32)), lf

    close(unit)

  end subroutine DVtkArrayBin
  
end module vtkarray


