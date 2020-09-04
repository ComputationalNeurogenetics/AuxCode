s.data<-readRDS("scATAC/E14 vMB1_phase3.Rds")

# Run fix for having chrUn in peaks data (corrected above but not run again)
DefaultAssay(s.data) <- 'peaks'
feat.keep <- bins[grep(pattern="chrUn.*",x=rownames(s.data), invert=T)]
s.data<-subset(s.data, features=feat.keep)
