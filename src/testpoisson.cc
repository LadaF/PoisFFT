#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "poisfft.h"

int main(){
  const int ntimes = 2;
  const int ns[3] = {11, 12, 13};
  
  double ds[3];
  
  const int BCs[6] = {PoisFFT::Periodic, PoisFFT::Periodic,
                      PoisFFT::Periodic, PoisFFT::Periodic,
                      PoisFFT::Periodic, PoisFFT::Periodic};
  
  int i;
  
  for (i = 0; i<3; i++){
   ds[i] = 1.0 / ns[i];
  }
  
  double *arr = (double *) malloc(ns[0]*ns[1]*ns[2] * sizeof(double));
  double *RHS = (double *) malloc(ns[0]*ns[1]*ns[2] * sizeof(double));
  
  memset(arr, 0, ns[0]*ns[1]*ns[2] * sizeof(double));
  memset(RHS, 0, ns[0]*ns[1]*ns[2] * sizeof(double));
  
  PoisFFT::Solver<3, double> S(ns, ds, BCs);
  
  for (i = 0; i<ntimes; i++){
    S.execute(arr, RHS);
  }

  free(arr);
  arr = NULL;
  free(RHS);
  RHS = NULL;
}
