    type(solver), pointer :: f_D
    type(c_ptr), intent(out) :: D
    integer(c_int) :: nxyz(dims), BCs(2*dims)
    real(rp) :: dxyz(dims)
    integer(c_int), optional :: gnxyz(dims), offs(dims)
    type(c_ptr), value :: mpi_comm
    integer(c_int), value :: nthreads
    integer :: f_comm
    integer, parameter :: idxs(6) = [2,1,4,3,6,5]
    
#ifdef MPI
    f_comm = MPI_Comm_c2f(mpi_comm)
#endif

    allocate(f_D)
    
    if (nthreads <1) nthreads = 1

    if (present(gnxyz).and.present(offs)) then
      f_D = solver(int(nxyz(dims:1:-1)), &
                   dxyz(dims:1:-1), &
                   int(BCs(idxs(2*dims:1:-1))), &
                   int(gnxyz(dims:1:-1)), &
                   int(offs(dims:1:-1)), &
                   f_comm, &
                   int(nthreads))
    else
      f_D = solver(int(nxyz(dims:1:-1)), &
                   dxyz(dims:1:-1), &
                   int(BCs(idxs(2*dims:1:-1))), &
                   nthreads=int(nthreads))
    end if

    D = c_loc(f_D)
