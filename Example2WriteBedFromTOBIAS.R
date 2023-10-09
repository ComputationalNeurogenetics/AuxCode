library(qs)
library(valr)
library(GenomicRanges)
library(devtools)
source_url("https://github.com/ComputationalNeurogenetics/AuxCode/raw/master/AuxFunctions.R")

# Read TOBIAS data as GRange
rV2.groups.tobias.gr <-qread("analysis/rV2.groups.tobias.gr.qs")

# Writes TOBIAS bound location of TF in group into bed file
ConstructBed_TobiasGr(gr=rV2.groups.tobias.gr, group="CO1",TF="Tal1", file="test.bed")