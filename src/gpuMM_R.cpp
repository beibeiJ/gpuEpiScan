//#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <R_ext/Rdynload.h>
#include "gpuMatrixMultiply.h"



extern "C"
SEXP gpuMM_R(SEXP A,
             SEXP B,
             SEXP m,
             SEXP k,
             SEXP n,
             SEXP idx){
  SEXP C;
  PROTECT(C = allocVector(REALSXP, asInteger(m)*asInteger(n)));
  gpumm(REAL(A), REAL(B), INTEGER(m), INTEGER(k), INTEGER(n), REAL(C), INTEGER(idx));
  UNPROTECT(1);
  return C;
}


void R_init_gpuEpiScan(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
