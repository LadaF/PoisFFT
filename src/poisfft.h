#ifdef __cplusplus
extern "C" {
#endif
const int POISFFT_PERIODIC = 0;
const int POISFFT_DIRICHLET = 1;
const int POISFFT_NEUMANN = 2;
const int POISFFT_DIRICHLET_STAG = 3;
const int POISFFT_NEUMANN_STAG = 4;
const int POISFFT_SPECTRAL = 0;
const int POISFFT_FINITE_DIFFERENCE_2 = 2;
 
typedef struct{int dims; void *D;} poisfft_solver;

typedef struct{int dims; void *D;} poisfft_solver_f;

poisfft_solver poisfft_solver_new(const int dims,
                         const int *nxyz, const double *Lxyz,
                         const int *BCs, const int approximation, 
                         const int *gnxyz, const int *offs,
                         const void *mpi_comm, int nthreads);
void poisfft_solver_execute( poisfft_solver S, double *Phi, const double *RHS,
                             const int *ngPhi, const int *ngRHS);
void poisfft_solver_finalize( poisfft_solver *S);

poisfft_solver_f poisfft_solver_f_new(const int dims,
                         const int *nxyz, const float *Lxyz,
                         const int *BCs, const int approximation,
                         const int *gnxyz, const int *offs,
                         const void *mpi_comm, int nthreads);
void poisfft_solver_f_execute( poisfft_solver_f S, float *Phi, const float *RHS,
                             const int *ngPhi, const int *ngRHS);
void poisfft_solver_f_finalize( poisfft_solver_f *S);

#ifdef __cplusplus
}

namespace PoisFFT{
  const int PERIODIC = POISFFT_PERIODIC;
  const int DIRICHLET = POISFFT_DIRICHLET;
  const int NEUMANN = POISFFT_NEUMANN;
  const int DIRICHLET_STAG = POISFFT_DIRICHLET_STAG;
  const int NEUMANN_STAG = POISFFT_NEUMANN_STAG;
  const int SPECTRAL = POISFFT_SPECTRAL;
  const int FINITE_DIFFERENCE_2 = POISFFT_FINITE_DIFFERENCE_2;

  template <unsigned int dims, typename real>  class Solver{
    private:
    public:
      Solver(const int *nxyz, const real *Lxyz,
                     const int *BCs, const int *gnxyz=0, const int *offs=0,
                     const void *mpi_comm=0, int nthreads=1);
      ~Solver();
      void execute(real *Phi, const real *RHS,
                   const int *ngPhi=0, const int *ngRHS=0);
  };
  
  template <unsigned int dims>  class Solver<dims, double>{
    private:
      poisfft_solver c_solver;
    public:
      Solver(const int *nxyz, const double *Lxyz,
                     const int *BCs, const int approximation=0,
                     const int *gnxyz=0, const int *offs=0,
                     const void *mpi_comm=0, int nthreads=1){
        c_solver = poisfft_solver_new(dims, nxyz, Lxyz, BCs, approximation, 
                                      gnxyz, offs, mpi_comm, nthreads);
      };
      ~Solver(){
        poisfft_solver_finalize(&c_solver);
      };
      void execute(double *Phi, const double *RHS,
                   const int *ngPhi=0, const int *ngRHS=0){
        poisfft_solver_execute(c_solver, Phi, RHS, ngPhi, ngRHS);
      };
  };
  
  template <unsigned int dims>  class Solver<dims, float>{
    private:
      poisfft_solver_f c_solver;
    public:
      Solver(const int *nxyz, const float *Lxyz,
                     const int *BCs, const int approximation=0,
                     const int *gnxyz=0, const int *offs=0,
                     const void *mpi_comm=0, int nthreads=1){
        c_solver = poisfft_solver_f_new(dims, nxyz, Lxyz, BCs, approximation,
                                        gnxyz, offs, mpi_comm, nthreads);
      };
      ~Solver(){
        poisfft_solver_f_finalize(&c_solver);
      };
      void execute(float *Phi, const float *RHS,
                   const int *ngPhi=0, const int *ngRHS=0){
        poisfft_solver_f_execute(c_solver, Phi, RHS, ngPhi, ngRHS);
      };
  };

}
#endif
