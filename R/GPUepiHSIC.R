
##' @description Calculate the significance of epistasis according to the definition of HSIC, conduct Z test for HSIC values and
##' choose variant pairs with the significance below the given threshold for output.
##' @title Calculate epistasis on GPU using HSIC with two genotype inputs
##' @export
##' @param geno1 is the first normalized genotype data. It can be a matrix or a dataframe, or a big.matrix object from \pkg{bigmemory}.
##' The columns contain the information of variables and the rows contain the information of samples.
##' @param geno2 is the second normalized genotype data. It can be a matrix or a dataframe, or a big.matrix object from \pkg{bigmemory}.
##' The columns contain the information of variables and the rows contain the information of samples.
##' @param pheno is a vector containing the normalized phenotype information.
##' @param zpthres is the significance threshold for cut-off output of the variant pairs.
##' @param chunk is the number of variants in each chunk.
##' @param outfile is the prefix of out filename.
##' @param suffix is the suffix of out filename.
##' @param gpuidx is the integer index of GPU(device) starting from 0. The default is 0.
##' @param ... not used
##' @import stats
##' @return null
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' # simulate some data
#' set.seed(123)
#' geno1 <- matrix(sample(0:2, size = 1000, replace = TRUE, prob = c(0.5, 0.3, 0.2)), ncol = 10)
#' geno2 <- matrix(sample(0:2, size = 2000, replace = TRUE, prob = c(0.4, 0.3, 0.3)), ncol = 20)
#' dimnames(geno1) <- list(row = paste0("IND", 1:nrow(geno1)), col = paste0("rs", 1:ncol(geno1)))
#' dimnames(geno2) <- list(row = paste0("IND", 1:nrow(geno2)), col = paste0("exm", 1:ncol(geno2)))
#' p2 <- rnorm(100, mean = 5, sd = 10)
#'
#' # normalized data
#' geno1 <- scale(geno1)
#' geno2 <- scale(geno2)
#' p2 <- as.vector(unlist(scale(p2)))
#'
##' # two genotypes with quantitative phenotype
#' GPUepiHSIC2genos(geno1 = geno1,
#' geno2 = geno2,
#' pheno = p2,
#' outfile = "episcan_2geno_quant",
#' suffix = ".txt",
#' zpthres = 0.1,
#' chunk = 10)
#'
#' # take a look at the result
#' res <- read.table("episcan_2geno_quant.txt", header = TRUE, stringsAsFactors = FALSE)
#' head(res)
GPUepiHSIC2genos <- function(geno1 = NULL,
                             geno2 = NULL,
                             pheno = NULL,
                             chunk = 1000,
                             zpthres = 10e-06,
                             outfile = "NONE",
                             suffix = ".txt",
                             gpuidx = 0,
                             ...){
  zthres <- abs(qnorm(zpthres/2))
  ## output head
  # check whether output file exisit or not; if yes, delete
  OUT <- paste(outfile, suffix, sep = "")
  if(file.exists(OUT)) file.remove(OUT)
  cat(paste("SNP1", "SNP2", "Zscore", "ZP",
            sep = " "),
      "\n",
      file = OUT,
      append = TRUE)
  nSNP1 <- ncol(geno1)
  nsplits1 <- ceiling(nSNP1/chunk)

  nSNP2 <- ncol(geno2)
  nsplits2 <- ceiling(nSNP2/chunk)


  gc()
  # first loop for all the ones which is interger times than chunk
  for ( i in 1:nsplits1)
  {
    print(paste(i, "chunk loop:", date()))
    for (j in 1:nsplits2)
    {
      HSIC.Zmatrix <- GPUepiHSIC(A = geno1[, ithChunk(i, nSNP1, chunk), drop = FALSE],
                                 B = geno2[, ithChunk(j, nSNP2, chunk), drop = FALSE],
                                 P = pheno,
                                 GPUidx = gpuidx)
      index <- which(abs(HSIC.Zmatrix) >= zthres, arr.ind = TRUE)
      WriteSnpPairs(Zmatrix = HSIC.Zmatrix, indexArr = index,
                    outfile = OUT)
      rm(list = c("HSIC.Zmatrix", "index"))
      gc()
    }
  }


  # finish
  print("GPUepiHSIC calculation is over!")
  print(date())
}











##' @description Calculate HSIC values on GPU
##' @title Calculate HSIC values on GPU
##' @export
##' @param A is one matrix.
##' @param B is one matrix.
##' @param P is "phenoype", a vector.
##' @param s is the standard deviation of the epiHSIC distribution
##' @param GPUidx is the integer index of GPU(device) starting from 0. Default is 0.
##' @param ... not used.
##' @return a matrix
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
GPUepiHSIC <- function(A = NULL,
                       B = NULL,
                       P = NULL,
                       s = NULL,
                       GPUidx = 0,
                       ...)
{
  if(is.null(s)) s <- sqrt(length(P))
  H.Zmatrix <- gpuMM(A = t(A),
                     B = P * B,
                     m = ncol(A),
                     k = nrow(A),
                     n = ncol(B),
                     idx = GPUidx) / s
  rownames(H.Zmatrix) <- colnames(A)
  colnames(H.Zmatrix) <- colnames(B)
  return(H.Zmatrix)
}










##' @description Calculate the significance of epistasis according the definition of HSIC, conduct Z test for HSIC values and
##' choose variant pairs with the significance below the given threshold for output.
##' @title Calculate epistasis on GPU using HSIC with one genotype input
##' @export
##' @param geno is the normalized genotype data. It can be a matrix or a dataframe, or a big.matrix object from \pkg{bigmemory}.
##' The columns contain the information of variables and the rows contain the information of samples.
##' @param pheno is a vector containing the normalized phenotype information.
##' @param zpthres is is the significance threshold to select variant pairs for output. Default is 1e-6.
##' @param chunk is the number of variants in each chunk.
##' @param outfile is the prefix of out filename.
##' @param suffix is the suffix of out filename.
##' @param gpuidx is the integer index of GPU(device) starting from 0. The default is 0.
##' @param ... not used.
##' @import stats
##' @return null
##' @author Beibei Jiang \email{beibei_jiang@@psych.mpg.de}
##' @examples
##' # simulate some data
#' set.seed(123)
#' geno1 <- matrix(sample(0:2, size = 1000, replace = TRUE, prob = c(0.5, 0.3, 0.2)), ncol = 10)
#' dimnames(geno1) <- list(row = paste0("IND", 1:nrow(geno1)), col = paste0("rs", 1:ncol(geno1)))
#' p2 <- rnorm(100, mean = 5, sd = 10)
#'
#' # normalized data
#' geno1 <- scale(geno1)
#' p2 <- as.vector(unlist(scale(p2)))
#'
#' # one genotypes with quantitative phenotype
#' GPUepiHSIC1geno(geno = geno1,
#' pheno = p2,
#' outfile = "episcan_1geno_quant",
#' suffix = ".txt",
#' zpthres = 0.1,
#' chunk = 10)
#'
#' # take a look at the result
#' res <- read.table("episcan_1geno_quant.txt", header = TRUE, stringsAsFactors = FALSE)
#' head(res)
GPUepiHSIC1geno <- function(geno = NULL,
                            pheno,
                            chunk = 1000,
                            zpthres = 10e-06,
                            outfile = "NONE",
                            suffix = ".txt",
                            gpuidx = 0,
                            ...)
{
  zthres <- abs(qnorm(zpthres/2))
  ## output head
  # check whether output file exisit or not; if yes, delete
  OUT <- paste(outfile, suffix, sep = "")
  if(file.exists(OUT)) file.remove(OUT)
  cat(paste("SNP1", "SNP2", "Zscore", "ZP",
            sep = " "),
      "\n",
      file = OUT,
      append = TRUE)

  nSNP <- ncol(geno)
  nsplits <- ceiling(nSNP/chunk)
  gc()

  # first loop for all the ones which is interger times than chunk
  for ( i in 1:nsplits)
  {
    print(paste(i, "chunk loop:", date()))
    for (j in i:nsplits)
    {
      HSIC.Zmatrix <- GPUepiHSIC(A = as.matrix(geno[, ithChunk(i, nSNP, chunk), drop = FALSE]),
                                 B = as.matrix(geno[, ithChunk(j, nSNP, chunk), drop = FALSE]),
                                 P = pheno,
                                 GPUidx = gpuidx)
      index <- which(abs(HSIC.Zmatrix) >= zthres, arr.ind = TRUE)
      ifelse(i == j,
             WriteSnpPairs_sym,
             WriteSnpPairs)(Zmatrix = HSIC.Zmatrix, indexArr = index,
                            outfile = OUT)
      rm(list = c("HSIC.Zmatrix", "index"))
      gc()
    }
  }



  print("GPUepiHSIC calculation is over!")
  print(date())
}


