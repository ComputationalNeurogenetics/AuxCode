library(GenomicRanges)
library(devtools)
library(parallel)
library(tidyverse)
source_url("https://github.com/ComputationalNeurogenetics/AuxCode/raw/master/AuxFunctions.R")

# Read TOBIAS results from snakemake output
rV2.groups.tobias <- get_BINDetect_snakemake_results_v2("/Volumes/LaCie/E12rV2_groups_200923/TFBS/", parallel = TRUE, mc.cores = 4)
