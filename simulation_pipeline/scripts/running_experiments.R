if(!require(BiocManager)) install.packages("BiocManager", repo = 'http://cran.us.r-project.org')
if(!require(polyester)) BiocManager::install("polyester")

library(polyester)
library(Biostrings)

fasta_file <- paste('', snakemake@input[["transcriptFasta"]], sep = '')

df <- read.csv(paste('', snakemake@input[["transcriptReadCounts"]], sep = ''), header = FALSE)

read_counts <- matrix(rep(df[[1]], each = snakemake@params[["nreplicates"]]), ncol = snakemake@params[["nreplicates"]], byrow = TRUE)


simulate_experiment_countmat(fasta_file, readmat=read_counts, outdir=paste('', snakemake@output[["outDir"]], sep = ''))
