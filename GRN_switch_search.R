GRN_switch_search <- function(GRN.data, seurat.object){
  require(Seurat)
  require(Cicero)
  require(tidyverse)
  require(presto)
  require(Signac)
  require(BSgenome.Mmusculus.UCSC.mm10)
  
  # Entire search space in terms of genes
  G.a <- unique(GRN.data$markers_rna$RNA.feature)
  
  # Prefilter peaks to G to links
  Pe.to.G <- GRN.data$peaks.to.genes[GRN.data$peaks.to.genes$score >= .35 & GRN.data$peaks.to.genes$pvalue <= 0.05,]
  Pe.to.G.posneg <- GRN.data$peaks.to.genes[(GRN.data$peaks.to.genes$score >= .35 | GRN.data$peaks.to.genes$score <= -0.35 )& GRN.data$peaks.to.genes$pvalue <= 0.05,]
  
  # Output formatting

  output <- lapply(G.a, function(G){
    G.output<-list()
    # Set the gene type either G.tf or G.non.tf
    if (G %in% GRN.data$markers_motifs$gene) {G.output$type <- "G.tf"} else {G.output$type <- "G.non.tf"}
    
    # If G.tf then gather G.tf expression to G.tf.z-score association for that G.tf, these are potential downstream targets
    if (G.output$type == "G.tf" & G %in% GRN.data$TF_RNA_TF_z_pairs$gene){G.output$downstream_AUC <- filter(GRN.data$TF_RNA_TF_z_pairs, gene == G)}
    
    # Looking for upstream of G what peaks are associated with expression
    if (G %in% Pe.to.G$gene) {
      G.output$Pe.to.G <- Pe.to.G[Pe.to.G$gene == G,]
      Pe.anno <- ClosestFeature(object = seurat.object, regions = G.output$Pe.to.G)
      mcols(G.output$Pe.to.G)<-cbind(mcols(G.output$Pe.to.G),Pe.anno)
      }
    
    # Then back downstream analysis
    if (G.output$type == "G.tf" & nrow(G.output$downstream_AUC)>0){
      # We need search through motif matrix that which genes are downstream_AUC targets
      G.tf.motif <- gsub(x = pull(G.output$downstream_AUC,motif.feature), replacement = "_", pattern = "-")
      G.output$downstream_targets <- GRN.data$motifs@positions[G.tf.motif][[1]]
      # Let's annotate downstream targets with the logic of closest known feature
      Downstream.anno <- ClosestFeature(object = seurat.object, regions = G.output$downstream_targets)
      mcols(G.output$downstream_targets)<-cbind(mcols(G.output$downstream_targets),Downstream.anno)
      # Here I need to look RNA marker data from GRN.data$markers_rna, and compare it to the downstream gene list. 
      # How to make some statistical analysis of it?
      # At least question could be turned around, which of the rna marker genes are "explained" via link to any GRN switch
      # NOTE, this list of potential downstream genes can be further filtered with the same peaks to genes calculation, and possible even further by with footprinting
      
      # Then as filter step 1, we can filter potential downstream targets based on the fact that if this G.tf is associated in its expression control there needs to be 
      # some correlation (pos or neg) between the peak where the motif is and the gene. Same Peaks to Genes to association data is used here.
 
      G.downstream.filt1 <- G.output$downstream_targets[G.output$downstream_targets$gene_name %in% Pe.to.G.posneg$gene,]
      # Then filtering step 2 could be in the form of Footprinting as we should see the G.tf in question having some footprint in downstream gene promoter if it is actively regulating its expression
      Footprint(object = seurat.object, genome = BSgenome.Mmusculus.UCSC.mm10, regions = G.downstream.filt1, key="Dlx2foot")
      p2 <- PlotFootprint(tmp, features = c("Dlx2foot"), idents = c("P3 GABAergic neurons", "P1 GABA"))
      
      # Here problem with the footprinting is interpretation of the results, and that it is still aggregate over all siten given in the regions argument
      # Filtering step 3 could be filtering through tissue specific marker list
      
      }
    
    
    # Then we need CRE information in more usable format, perhaps no ccans yet
    
    # We need footprinting to P
    
    return(G.output)
  })
  
  names(output) <- G.a
  
}

format_GRN_search_data<-function(seurat.object, ident.of.interest, cell.pseudotime = NULL, padj.cutoff = 1e-2, CRE.score.cutoff = .5, txdb = TxDb.Mmusculus.UCSC.mm10.ensGene){
  # Assuming now that seurat.object contains assays peaks, RNA, chromvar, also links and motif.match for footprinting
  
  require(Seurat)
  require(Cicero)
  require(tidyverse)
  require(presto)
  require(Signac)
  
  # Initialize output object
  output <- list()
  
  # Write here smoothing with pseudo time if pseudotime has been given
  
  #---------------------------
  
  # Using presto based WX test to find TFs which are both expressed and which motif are accessible
  print("Running presto::wilcoxauc for RNA assay")
  markers_rna <- presto:::wilcoxauc.Seurat(X = seurat.object, group_by = 'predicted.id', assay = 'data', seurat_assay = 'RNA')
  print("Running presto::wilcoxauc for chromvar assay")
  markers_motifs <- presto:::wilcoxauc.Seurat(X = seurat.object, group_by = 'predicted.id', assay = 'data', seurat_assay = 'chromvar')
  
  # Formatting resulting objects
  motif.names <- markers_motifs$feature
  colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
  colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
  markers_rna$gene <- markers_rna$RNA.feature
  # ConverMotifID does not work with 10x motifs due to the lacking mapping between motif id and name (both catenated to string)
  #markers_motifs$gene <- ConvertMotifID(seurat.object, id = motif.names)
  # Custom code to extract gene names, remove (var) infor, and change to capital first letter format
  markers_motifs$gene <-  str_to_title(sapply(strsplit(motif.names,"-"), function(x){x[1]})) %>% gsub(pattern = "\\(.*\\)", replacement = "")
  
  # Call for separate function (from Satija lab pages) for detection of interesting TF and TF-motif pairs
  print("Find TF_RNA and TF_z-scorea associations")
  TF_RNA_TF_z_pairs <- topTFs(celltype = ident.of.interest, padj.cutoff = padj.cutoff)
  
  output$markers_rna <- markers_rna
  output$markers_motifs <- markers_motifs
  output$TF_RNA_TF_z_pairs <- tibble(TF_RNA_TF_z_pairs)
  
  #----------------------------
  
  # Analyze through Cicero CRE links
  CRE.data <- Links(seurat.object)
  filt.CRE.data <- CRE.data[CRE.data$score >= CRE.score.cutoff]
  # Let's annotate CRE data, it is granges object thus let's annotate it as such
  
  txdb_str_chrom <- keepStandardChromosomes(txdb)
  mm10.promoters<- promoters(txdb_str_chrom, upstream = 2000, downstream = 400)
  
  # Finding which CREs are in promoter area and which not
  promoter.CRE <- filt.CRE.data[filt.CRE.data %over% mm10.promoters]
  non.promoter.CRE <- filt.CRE.data[filt.CRE.data %outside% mm10.promoters]
  peaks.to.genes <- LinkPeaks(s.data, peak.assay = "peaks", expression.assay = "RNA")

  output$promoter.CRE <- promoter.CRE
  output$non.promoter.CRE <- non.promoter.CRE
  output$peaks.to.genes <- Links(peaks.to.genes)

  output$motifs <- Motifs(seurat.object)
  
  return(output)
}
 
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}