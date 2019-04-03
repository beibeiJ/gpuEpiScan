# gpuEpiScan
Bug report git repository for gpuEpiScan


## Introduction

**gpuEpiScan** provides efficient methods, which were implemented based on GPU, to scan for pairwise epistasis in both case-control study and quantitative studies. It is suitable for genome-wide interaction studies (GWIS) by splitting the computation into manageable chunks. The epistasis methods used by **gpuEpiScan** are adjusted from two published papers (see References below).
  
  Unfortunately, the current **gpuEpiScan** version only support Linux system. We are trying to figure out how to make it work with Windows system. 
  
  ## Installation
  To install **gpuEpiScan** sucessfully, `CUDA` and `nvcc` must be installed. If your `CUDA` path is not `/usr/local/cuda` and your `nvcc` is not global command, a specific `CUDA` path needs to be given by the user.
  
```{r, eval = F, error = F, results='hide'}
# normal install, the default will take CUDA path with /usr/local/cuda 
# or the detected path with shell command 'which nvcc'
devtools::install_git('https://github.com/beibeiJ/gpuEpiScan.git')

# if the above installation doesn't work, `CUDA` path needs to be specified
devtools::install_git('https://github.com/beibeiJ/gpuEpiScan.git', 
                      build_opts = '--configure-args=YOUR_CUDA_PATH')
```

An alternative way is to download the source code of **gpuEpiScan** then install via

```{r, eval = F, error = F, results='hide'}
# x.y.z need to be replaced by the real version, e.g., 0.0.1
R CMD INSTALL gpuEpiScan_x.y.z.tar.gz --configure-args=YOUR_CUDA_PATH
```


## Sample implementation
```{r, echo = T, results='hide'}
# load package
library(gpuEpiScan)
```
First, we generate a small genotype dataset (`geno`) with sample size of 100 subjects and 1000 variables (e.g., SNPs) as well as a case-control phenotype (`p`).  
```{r, echo = T, results='markup'}
set.seed(123)
geno <- matrix(sample(0:2, size = 100000, 
                      replace = TRUE, prob = c(0.5, 0.3, 0.2)), ncol = 1000)
dimnames(geno) <- list(row = paste0("IND", 1:nrow(geno)), 
                       col = paste0("rs", 1:ncol(geno)))
p <- c(rep(0, 60), rep(1, 40))
geno[1:5, 1:8]
```
To  use **gpuEpiScan**, simply start with the main function `gpuEpiScan`. There are three mandatory parameters to be set by the user: genotype data, phenotype data and phenotype category ("case-control" or "quantitative"). Since the data simulated above is not normalized yet, we need to set `scale = TRUE`. By passing an integer number to parameter `chunksize`, the genotype data will be split into several chunks of that size during the calculation. For the example above (using `geno`), `chunksize = 200` means each chunk has 20 variables(variants) and the total number of chunks is 5. Moreover, in most cases, the result of epistasis analysis is huge due to the large number of the variable (variants) combinations. To reduce the size of the result file, setting a threshold of the statistical test (`zpthres`) to have an output cut-off is a practical option. By seting `gpuidx = 0`, the GPU calculation will run on the GPU device whose index number is equal to 0.
```{r, echo = T, results='markup'}
gpuEpiScan(geno1 = geno, 
           pheno = p, 
           phetype = "case-control",
           outfile = "episcan_1geno_cc", 
           suffix = ".txt", 
           zpthres = 0.9, 
           chunksize = 200, 
           scale = TRUE,
           gpuidx = 0)
```
The result of `gpuEpiScan` is stored in the specified file ("episcan_1geno_cc.txt"). Let's take a look:
```{r, echo = T, results='markup'}
result <- read.table("episcan_1geno_cc.txt",
                     header = TRUE,
                     stringsAsFactors = FALSE)
head(result)
```
In a genome-wide level epistasis study, it is usual to have millions of variables (variants). Analyzing such big data is super time-consuming. The common appoach is to parallelize the task and run the subtasks with High Performance Computing (HPC) techniques, e.g., on a cluster. By splitting genotype data per chromosome, the huge epistasis analysis task can be divided into relatively small tasks. If only 22 chromosomes exist in the initial task, there are 253 ((1+22)\*22/2) subtasks after splitting and considering all the combinations of the chromosomes. `gpuEpiScan` supports two inputs of genotype data by simply passing information to "geno1" and "geno2".  
```{r, echo = T, results='markup'}
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
```


## References
- Kam-Thong, T., D. Czamara, K. Tsuda, K. Borgwardt, C. M. Lewis, A. Erhardt-Lehmann, B. Hemmer, et al. (2011). doi:10.1038/ejhg.2010.196.

- Kam-Thong, T., B. Pütz, N. Karbalai, B. Müller-Myhsok, and K. Borgwardt. (2011). doi:10.1093/bioinformatics/btr218
