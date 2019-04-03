## ----setup, include=FALSE------------------------------------------------
  knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
  options(stringsAsFactors = FALSE)

## ---- eval = F, error = F, results='hide'--------------------------------
#  # normal install, the default will take CUDA path with /usr/local/cuda
#  # or the detected path with shell command 'which nvcc'
#  devtools::install_git('https://github.com/beibeiJ/gpuEpiScan.git')
#  
#  # if the above installation doesn't work, `CUDA` path needs to be specified
#  devtools::install_git('https://github.com/beibeiJ/gpuEpiScan.git',
#                        build_opts = '--configure-args=YOUR_CUDA_PATH')

## ---- eval = F, error = F, results='hide'--------------------------------
#  # x.y.z need to be replaced by the real version, e.g., 0.0.1
#  R CMD INSTALL gpuEpiScan_x.y.z.tar.gz --configure-args=YOUR_CUDA_PATH

## ---- echo = T, results='hide'-------------------------------------------
# load package
library(gpuEpiScan)

## ---- echo = T, results='markup'-----------------------------------------
set.seed(123)
geno <- matrix(sample(0:2, size = 100000, 
                      replace = TRUE, prob = c(0.5, 0.3, 0.2)), ncol = 1000)
dimnames(geno) <- list(row = paste0("IND", 1:nrow(geno)), 
                       col = paste0("rs", 1:ncol(geno)))
p <- c(rep(0, 60), rep(1, 40))
geno[1:5, 1:8]

## ---- echo = T, results='markup'-----------------------------------------
gpuEpiScan(geno1 = geno, 
           pheno = p, 
           phetype = "case-control",
           outfile = "episcan_1geno_cc", 
           suffix = ".txt", 
           zpthres = 0.9, 
           chunksize = 200, 
           scale = TRUE,
           gpuidx = 0)

## ---- echo = T, results='markup'-----------------------------------------
result <- read.table("episcan_1geno_cc.txt",
                     header = TRUE,
                     stringsAsFactors = FALSE)
head(result)

## ---- echo = T, results='markup'-----------------------------------------
# simulate data
geno1 <- matrix(sample(0:2, size = 100000, 
                       replace = TRUE, prob = c(0.5, 0.3, 0.2)), ncol = 1000)
geno2 <- matrix(sample(0:2, size = 200000, 
                       replace = TRUE, prob = c(0.4, 0.3, 0.3)), ncol = 2000)
dimnames(geno1) <- list(row = paste0("IND", 1:nrow(geno1)), 
                        col = paste0("rs", 1:ncol(geno1)))
dimnames(geno2) <- list(row = paste0("IND", 1:nrow(geno2)), 
                        col = paste0("exm", 1:ncol(geno2)))
p <- rnorm(100)

# scan epistasis
gpuEpiScan(geno1 = geno1,
           geno2 = geno2,
           pheno = p, 
           phetype = "quantitative",
           outfile = "episcan_2geno_quant", 
           suffix = ".txt", 
           zpthres = 0.6, 
           chunksize = 200, 
           scale = TRUE)

