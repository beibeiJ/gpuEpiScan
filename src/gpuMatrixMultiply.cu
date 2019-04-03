/* The example is adapted from
https://github.com/sol-prog/cuda_cublas_curand_thrust */


#include <stdio.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "gpuMatrixMultiply.h"


/* Multiply the matrix A and matrix B on GPU and save the result into C */
/* C(m, n) = alpha * A(m, k) %*% B(k, n) + beta * C(m, n) */
/* If it is the matrix product of just A and B, the actual calculation would be
alpha = 1 and beta = 0, so C(m, n) = 1 * A(m, k) % *% B(k, n) + 0 * C(m, n) */

/* *A is the pointer which points to a matrix, same for *B and *C */

void gpu_blas_mm(double *A, double *B, double *C, int m, int k, int n)
{
	int nra = m, nrb = k, nrc = m; /* nr means number of rows */
	const double alf = 1.0;
	const double bet = 0.0;
	const double *alpha = &alf;
	const double *beta = &bet;

	/* create a product for cuBLAS */
	cublasHandle_t product;
	cublasCreate(&product);

	/* do the actual multiplication */
	cublasDgemm(product, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, nra, B, nrb, beta, C, nrc);

	/* destroy the product */
	cublasDestroy(product);

}


/* gpumm: this function performs matrix multiplication on GPU and send the result back to CPU */
/* C(m, n) = alpha * A(m, k) %*% B(k, n) + beta * C(m, n) */
void gpumm(double *A, double *B, int *m, int *k, int *n, double *C, int *idx)
{
  /* set working GPU(device) */
  cudaSetDevice(*idx);

	/* allocate 3 arrays on GPU(device) */
	double *gpuA, *gpuB, *gpuC;
	cudaMalloc(&gpuA, *m * *k * sizeof(double)); /* *m and *k is the value while m and k are the address */
	cudaMalloc(&gpuB, *k * *n * sizeof(double));
	cudaMalloc(&gpuC, *m * *n * sizeof(double));

	/* Copy CPU data to GPU */
	cudaMemcpy(gpuA, A, *m * *k * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(gpuB, B, *k * *n * sizeof(double), cudaMemcpyHostToDevice);


	/* matrix multiplication on GPU */
	gpu_blas_mm(gpuA, gpuB, gpuC, *m, *k, *n);

	/* copy the result from device to host memory */
	cudaMemcpy(C, gpuC, *m * *n * sizeof(double), cudaMemcpyDeviceToHost);

	/* free GPU memory */
	cudaFree(gpuA);
	cudaFree(gpuB);
	cudaFree(gpuC);
}


