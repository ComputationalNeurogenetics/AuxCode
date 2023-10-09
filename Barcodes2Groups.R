# Code for generating barcode to cluster mapping for TOBIAS runs
library(tidyverse)
library(Seurat)
library(Signac)
library(qs)

rV2.data <- qread("./scATAC_data/nmm_rV2_subset_relabeled_110522.qs", nthreads = 6)
barcodes2groups <- tibble(barcode=colnames(rV2.data),cluster=rV2.data$labeling)
write.table(barcodes2groups, file = "./analysis/barcodes2groups.csv",sep="\t",quote = FALSE,col.names = FALSE, row.names = FALSE)

