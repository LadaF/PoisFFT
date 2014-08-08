#include <cmath>
#include <iostream>
#include "poisfft.h"
#include <stdlib.h>

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
    std::cout << "OK\n" << std::endl;
  }
  else {
    std::cout << "FAIL, residuum" << sum << std::endl;
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
  const int BCs[6] = {PoisFFT::DIRICHLET_STAG, PoisFFT::DIRICHLET_STAG, 
                      PoisFFT::DIRICHLET_STAG, PoisFFT::DIRICHLET_STAG,
                      PoisFFT::DIRICHLET_STAG, PoisFFT::DIRICHLET_STAG};
  
  int i;
  
  // set the grid, depends on the boundary conditions
  for (i = 0; i<3; i++){
    ds[i] = Ls[i] / ns[i];
  }
  
  // allocate the arrays contiguously, you can use any other class
  // from which you can get a pointer to contiguous buffer
  double *arr = new double[ns[0]*ns[1]*ns[2]];
  double *RHS = new double[ns[0]*ns[1]*ns[2]];
  
  // set the right-hand side
  init_rhs(ns, ds, Ls, RHS);
  
  // create solver object, 3 dimensions, double precision
  PoisFFT::Solver<3, double> S(ns, Ls, BCs);
  
  //run the solver, can be run many times for different right-hand sides
  S.execute(arr, RHS);

  // check corectness
  check_solution(ns, ds, Ls, arr);

  delete[] RHS;
  delete[] arr;
}
