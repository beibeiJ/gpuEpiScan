#include <stdio.h>
void gpumm(double *A, double *B, int *m, int *k, int *n, double *C, int *idx);
void gpu_blas_mm(double *A, double *B, double *C, int m, int k, int n);
