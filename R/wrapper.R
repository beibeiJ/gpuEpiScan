#' Calculate matrix product on GPU: calling cublas library
#' @param A a matrix or a data.frame.
#' @param B a matrix or a data.frame.
#' @param m number of rows in A.
#' @param k number of columns in A. It should also be the number of rows in B.
#' @param n number of columns in B
#' @param idx the index of GPU(device). The default is 0.
#' @return a product matrix which has m rows and n columns.
#' @useDynLib gpuEpiScan, .registration = TRUE
#' @export
#'
#' @examples
#' set.seed(123)
#' m <- 100
#' k <- 200
#' n <- 100
#' A <- matrix(rnorm(m*k), ncol = k)
#' B <- matrix(rnorm(k*n), ncol = n)
#' gpuC <- gpuMM(A, B, m, k, n)
#' cpuC <- A%*%B
#' gpuC[1:6, 1:5]
#' cpuC[1:6, 1:5]

gpuMM <- function(A, B, m, k, n, idx = 0)
{
  res <- .Call("gpuMM_R",
               A = as.double(A),
               B = as.double(B),
               m = as.integer(m),
               k = as.integer(k),
               n = as.integer(n),
               idx = as.integer(idx))

  return(matrix(res, nrow = nrow(A), ncol = ncol(B)))
}

