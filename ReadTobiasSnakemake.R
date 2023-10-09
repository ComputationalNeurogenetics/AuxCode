library(GenomicRanges)
library(devtools)
library(parallel)
library(tidyverse)
library(qs)
source_url("https://github.com/ComputationalNeurogenetics/AuxCode/raw/master/AuxFunctions.R")

# Read TOBIAS results from snakemake output
rV2.groups.tobias.gr <- get_BINDetect_snakemake_results_gr("/Volumes/LaCie/E12rV2_groups_200923/TFBS/", parallel = TRUE, mc.cores = 4)

qsave(rV2.groups.tobias.gr, file = "analysis/rV2.groups.tobias.gr.qs", nthreads = 6)
