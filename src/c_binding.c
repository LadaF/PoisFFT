#include <stdio.h>


#define SOLVER(dims, suff) \
        typedef void *poisfft_solver ## dims ## d ## suff;

#define POIS_NEW(dims, suff, real) \
        extern void poisfft_solver ## dims ## d ## suff ## _new( \
                       poisfft_solver ## dims ## d ## suff *D,  \
                       const int *nxyz, const real *Lxyz, const int *BCs, \
                       const int *approximation, \
                       const int *gnxyz, const int *offs, const void *mpi_comm, \
                       int nthreads);
                          
#define POIS_EXECUTE(dims, suff, real) \
        extern void poisfft_solver ## dims ## d ## suff ## _execute( \
                       poisfft_solver ## dims ## d ## suff D,  \
                       real *Phi, const real *RHS, \
                       const int *ngPhi, const int *ngRHS);
                       
#define POIS_FINALIZE(dims, suff) \
        extern void poisfft_solver ## dims ## d ## suff ## _finalize( \
                       poisfft_solver ## dims ## d ## suff *D);

#define POIS_DIMS(dims, suff, real) \
  SOLVER(dims, suff)  \
  POIS_NEW(dims, suff, real) \
  POIS_EXECUTE(dims, suff, real) \
  POIS_FINALIZE(dims, suff)

#define POIS(suff, real) \
  POIS_DIMS(1, suff, real) \
  POIS_DIMS(2, suff, real) \
  POIS_DIMS(3, suff, real)

POIS(, double)

POIS(_f, float)

typedef struct{int dims; void *D;} poisfft_solver;

typedef struct{int dims; void *D;} poisfft_solver_f;

poisfft_solver poisfft_solver_new(const int dims,
                         const int *nxyz, const double *Lxyz, 
                         const int *BCs, const int approximation,
                         const int *gnxyz, const int *offs,
                         const void *mpi_comm, int nthreads){
 
    poisfft_solver S;
    
    S.dims = dims;
    
    S.D = NULL;

    switch (S.dims) {
      case 1:
        poisfft_solver1d_new((poisfft_solver1d *) &(S.D),
                             nxyz, Lxyz, BCs, &approximation, gnxyz, offs, mpi_comm, nthreads);
        break;
      case 2:
        poisfft_solver2d_new((poisfft_solver2d *) &(S.D), 
                             nxyz, Lxyz, BCs, &approximation, gnxyz, offs, mpi_comm, nthreads);
        break;
      case 3:
        poisfft_solver3d_new((poisfft_solver3d *) &(S.D),
                             nxyz, Lxyz, BCs, &approximation, gnxyz, offs, mpi_comm, nthreads);
        break;
      default:
        printf("Error, PoisFFT solver dimension must be 1, 2 or 3!\n");
        break;
    }
    
    return S;
}

void poisfft_solver_execute( poisfft_solver S, double *Phi, const double *RHS,
                             const int *ngPhi, const int *ngRHS){
 
    switch (S.dims) {
      case 1:
        poisfft_solver1d_execute((poisfft_solver1d) S.D, Phi, RHS, ngPhi, ngRHS);
        break;
      case 2:
        poisfft_solver2d_execute((poisfft_solver2d) S.D, Phi, RHS, ngPhi, ngRHS);
        break;
      case 3:
        poisfft_solver3d_execute((poisfft_solver3d) S.D, Phi, RHS, ngPhi, ngRHS);
        break;
      default:
        printf("Error, PoisFFT solver dimension must be 1, 2 or 3!\n");
        break;
    }
}

void poisfft_solver_finalize( poisfft_solver *S){
    switch (S->dims) {
      case 1:
        poisfft_solver1d_finalize((poisfft_solver1d *) &(S->D));
        break;
      case 2:
        poisfft_solver2d_finalize((poisfft_solver1d *) &(S->D));
        break;
      case 3:
        poisfft_solver3d_finalize((poisfft_solver1d *) &(S->D));
        break;
      default:
        printf("Error, PoisFFT solver dimension must be 1, 2 or 3!\n");
        break;
    }
}


void poisfft_solver_f_new(const int dims,
                         const int *nxyz, const float *Lxyz,
                         const int *BCs, const int *approximation,
                         const int *gnxyz, const int *offs,
                         const void *mpi_comm, int nthreads){
    poisfft_solver_f S;
 
    S.dims = dims;
    
    S.D = NULL;

    switch (S.dims) {
      case 1:
        poisfft_solver1d_f_new((poisfft_solver1d_f *) &(S.D),
                               nxyz, Lxyz, BCs, approximation,
                               gnxyz, offs, mpi_comm, nthreads);
        break;
      case 2:
        poisfft_solver2d_f_new((poisfft_solver2d_f *) &(S.D), 
                               nxyz, Lxyz, BCs, approximation,
                               gnxyz, offs, mpi_comm, nthreads);
        break;
      case 3:
        poisfft_solver3d_f_new((poisfft_solver3d_f *) &(S.D),
                               nxyz, Lxyz, BCs, approximation,
                               gnxyz, offs, mpi_comm, nthreads);
        break;
      default:
        printf("Error, PoisFFT solver dimension must be 1, 2 or 3!\n");
        break;
    }
}

void poisfft_solver_f_execute( poisfft_solver_f S, float *Phi, const float *RHS,
                               const int *ngPhi, const int *ngRHS){
 
    switch (S.dims) {
      case 1:
        poisfft_solver1d_f_execute((poisfft_solver1d_f) S.D, Phi, RHS, ngPhi, ngRHS);
        break;
      case 2:
        poisfft_solver2d_f_execute((poisfft_solver2d_f) S.D, Phi, RHS, ngPhi, ngRHS);
        break;
      case 3:
        poisfft_solver3d_f_execute((poisfft_solver3d_f) S.D, Phi, RHS, ngPhi, ngRHS);
        break;
      default:
        printf("Error, PoisFFT solver dimension must be 1, 2 or 3!\n");
        break;
    }
}

void poisfft_solver_f_finalize( poisfft_solver *S){
    switch (S->dims) {
      case 1:
        poisfft_solver1d_f_finalize((poisfft_solver1d_f *) &(S->D));
        break;
      case 2:
        poisfft_solver2d_f_finalize((poisfft_solver1d_f *) &(S->D));
        break;
      case 3:
        poisfft_solver3d_f_finalize((poisfft_solver1d_f *) &(S->D));
        break;
      default:
        printf("Error, PoisFFT solver dimension must be 1, 2 or 3!\n");
        break;
    }
}


