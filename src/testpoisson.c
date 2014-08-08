#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "poisfft.h"

const double pi = 3.14159265358979323846;

#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k

void init_rhs(const int ns[3], const double ds[3], const double Ls[3], double* a){
  int i,j,k;
  for (i=0;i<ns[0];i++){
    for (j=0;j<ns[1];j++){
      for (k=0;k<ns[2];k++){
        double x = ds[0]*(i+0.5);
        double y = ds[1]*(j+0.5);
        double z = ds[2]*(k+0.5);
        a[IND(i,j,k)] = sin(3*pi*x/Ls[0]) * sin(5*pi*y/Ls[1]) *sin(7*pi*z/Ls[2]);
      }
    }
  }
}

void check_solution(const int ns[3], const double ds[3], const double Ls[3], double* a){
  int i,j,k;
  double sum = 0;
  for (i=0;i<ns[0];i++){
    for (j=0;j<ns[1];j++){
      for (k=0;k<ns[2];k++){
        double x = ds[0]*(i+0.5);
        double y = ds[1]*(j+0.5);
        double z = ds[2]*(k+0.5);
        sum += pow(a[IND(i,j,k)] - 
                   (- 1.0 / ( pow(3*pi/Ls[0], 2.0) + 
                              pow(5*pi/Ls[1], 2.0) + 
                              pow(7*pi/Ls[2], 2.0) ) *
                      sin(3*pi*x/Ls[0]) * 
                      sin(5*pi*y/Ls[1]) *
                      sin(7*pi*z/Ls[2]))
                   , 2.0);
      }
    }
  }
  if (sum < 1e-10) {
    printf("OK\n");
  }
  else {
    printf("FAIL, residuum %f\n", sum);
  }
}

int main(){
  // domain dimensions
  const double Ls[3] = {1.1, 2.0, 2.9};

  // gridpoint numbers
  const int ns[3] = {71, 93, 101};
  
  // distances between gridpoints
  double ds[3];
  
  //boundary conditions
  const int BCs[6] = {POISFFT_DIRICHLET_STAG, POISFFT_DIRICHLET_STAG, 
                      POISFFT_DIRICHLET_STAG, POISFFT_DIRICHLET_STAG,
                      POISFFT_DIRICHLET_STAG, POISFFT_DIRICHLET_STAG};
  
  int i;
  
  // set the grid, depends on the boundary conditions
  for (i = 0; i<3; i++){
    ds[i] = Ls[i] / ns[i];
  }

  // allocate the arrays contiguously
  double *arr = (double *) malloc(ns[0]*ns[1]*ns[2] * sizeof(double));
  double *RHS = (double *) malloc(ns[0]*ns[1]*ns[2] * sizeof(double));

  // set the right-hand side
  init_rhs(ns, ds, Ls, RHS);

  // create solver object
  poisfft_solver S = poisfft_solver_new(3, ns, Ls, BCs, POISFFT_SPECTRAL,
                                        NULL, NULL, NULL, 0);

  //run the solver, can be run many times for different right-hand sides
  poisfft_solver_execute(S, arr, RHS, NULL, NULL);

  // free the memory of the solver
  poisfft_solver_finalize(&S);

  // check corectness
  check_solution(ns, ds, Ls, arr);

  free(arr);
  arr = NULL;
  free(RHS);
  RHS = NULL;
}