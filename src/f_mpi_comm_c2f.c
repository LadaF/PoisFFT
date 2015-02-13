#ifdef MPI

#include <mpi.h>

// This function is callable from Fortran. MPI_Comm_c2f itself may be just a macro.

MPI_Fint f_MPI_Comm_c2f(MPI_Comm *comm) {
  return MPI_Comm_c2f(*comm);
}

#endif

