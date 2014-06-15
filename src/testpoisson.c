#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "poisfft.h"

int main(){
  const int ntimes = 2;
  const int ns[3] = {11, 12, 13};
  
  double ds[3];
  
  const int BCs[6] = {PoisFFT_NeumannStag, PoisFFT_NeumannStag, 
                      PoisFFT_Periodic, PoisFFT_Periodic,
                      PoisFFT_Periodic, PoisFFT_Periodic};
  
  int i;
  
  for (i = 0; i<3; i++){
   ds[i] = 1.0 / ns[i];
  }
  
  double *arr = (double *) malloc(ns[0]*ns[1]*ns[2] * sizeof(double));
  double *RHS = (double *) malloc(ns[0]*ns[1]*ns[2] * sizeof(double));
  
  memset(arr, 0, ns[0]*ns[1]*ns[2] * sizeof(double));
  memset(RHS, 0, ns[0]*ns[1]*ns[2] * sizeof(double));
  
  poisfft_solver S = poisfft_solver_new(3, ns, ds, BCs, NULL, NULL, NULL, 0);
  
  for (i = 0; i<ntimes; i++){
    poisfft_solver_execute(S, arr, RHS, NULL, NULL);
  }
  
  poisfft_solver_finalize(&S);
  
  free(arr);
  arr = NULL;
  free(RHS);
  RHS = NULL;
}