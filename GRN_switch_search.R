is.RNA.marker <- function(G, GRN.data,RNA.padj.cutoff =  0.05, RNA.logFC = .5){
  markers <- filter(GRN.data$markers_rna, RNA.padj.cutoff <= RNA.padj.cutoff, RNA.logFC >= RNA.logFC) %>% pull(gene)
  is.marker <- G %in% markers
  return(is.marker)
}

is.G.tf <- function(G, tf.genes){
  if (G %in% tf.genes){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.Pe.to.Grna <- function(G, Pe.to.G.posneg, seurat.object){
  if (G %in% Pe.to.G.posneg$gene) {
    Pe.to.G <- Pe.to.G.posneg[Pe.to.G.posneg$gene == G,]
    DefaultAssay(seurat.object) <- 'peaks'
    Pe.anno <- ClosestFeature(object = seurat.object, regions = Pe.to.G)
    mcols(Pe.to.G)<-cbind(mcols(Pe.to.G), Pe.anno)
  } else {
    # If no Pe to G expression link found
    Pe.to.G <- NA
  }
  return(Pe.to.G)
}

G.tf.downstream.targets <- function(G,GRN.data, seurat.object){
    downstream_AUC <- filter(GRN.data$TF_RNA_TF_z_pairs, gene == G)
    # We need to search through motif matrix that which genes are downstream_AUC targets
    G.tf.motif <- gsub(x = pull(downstream_AUC,motif.feature), replacement = "_", pattern = "-")
    downstream_targets <- GRN.data$motifs@positions[G.tf.motif][[1]]
    # Let's annotate downstream targets with the logic of closest known feature
    Downstream.anno <- ClosestFeature(object = seurat.object, regions = downstream_targets)
    mcols(downstream_targets)<-cbind(mcols(downstream_targets),Downstream.anno)

    # Then as filter step 1, we can filter potential downstream targets based on the fact that if this G.tf is associated in its expression control there needs to be
    # some correlation (pos or neg) between the peak where the motif is and the gene. Same Peaks to Genes to association data is used here.

    # Let's see which genes having motif in question in the peaks they have correlation with (either pos or neg)
    G.downstream.filt1 <- downstream_targets[downstream_targets %over% Pe.to.G.posneg]
    # Following tests only based on gene names, above one should be more accurate, but not certain
    #G.downstream.filt1 <- G.output$downstream_targets[G.output$downstream_targets$gene_name %in% Pe.to.G.posneg$gene,]

    # Filtering step 2 could be filtering through tissue specific marker list
    # Acute TODO: is sorting out into what extent this should run as focused on single cluster vs. all
    G.downstream.filt2 <- G.downstream.filt1[G.downstream.filt1$gene_name %in% RNA.markers.ident]

    downstream_filtered <- tibble(gene_name = G.downstream.filt1$gene_name, RNA_marker = ifelse(G.downstream.filt1$gene_name %in% RNA.markers.ident, TRUE, FALSE), TF = ifelse(G.downstream.filt1$gene_name %in% GRN.data$markers_motifs$gene, TRUE, FALSE)) %>% distinct()
    return(downstream_filtered)
}

fetch.associated.CREs <- function(GRN.data, CRE.score.cutoff, enhancers, G.output){
  # CREs have already been sorted to promoter and non.promoter ones throughout the genome
  # Annotating these results based on if CRE happens to be in Pe linked to G(rna)
  #browser()
  # Lets first take promoter region coordinates for G, there could be multiple ones as some transcript could have different TSS(?)
  G.ensg.promoter <- GRN.data$promoter.coordinates[unlist(mcols(GRN.data$promoter.coordinates)$GENEID == G.ensg)]
  
  # As we are now looking data at gene level we can collapse all overlapping ranges
  G.ensg.promoter <- GenomicRanges::reduce(G.ensg.promoter)
  
  # Lets take all cicero recognized CRE elements overlapping with the promoter coordinate(s) of G
  G.promoter.CRE <- GRN.data$cicero_links[GRN.data$cicero_links %over% G.ensg.promoter & GRN.data$cicero_links$score >= CRE.score.cutoff]
  
  # Lets see to which ccan groups these CREs belong to
  G.ccan.groups <- unique(unlist(mcols(G.promoter.CRE)$group))
  
  # Lets search for all CREs in these groups
  G.non.promoter.CRE <- GRN.data$cicero_links[mcols(GRN.data$cicero_links)$group %in% G.ccan.groups & GRN.data$cicero_links$score >= CRE.score.cutoff]
  
  # Then lets further annotate CREs related to G with information whether any G.rna associated Pe overlaps with them
  # This is obviously valid only if there is Pe.to.G.rna links
  if (class(G.output$Pe.to.G) == "GRanges"){
    mcols(G.promoter.CRE, level="within")$G.rna.association <- G.promoter.CRE %over% G.output$Pe.to.G
  }
  
  # Then lets further annotate CREs based on whether they are in mouse enhancer atlas
  if (class(G.promoter.CRE) == "GRanges"){
    mcols(G.promoter.CRE, level="within")$known <- G.promoter.CRE %over% enhancers
  }
  
  if (class(G.non.promoter.CRE) == "GRanges"){
    mcols(G.non.promoter.CRE, level="within")$known <- G.non.promoter.CRE %over% enhancers
  }
  return(list(G.promoter.CRE=G.promoter.CRE, G.non.promoter.CRE=G.non.promoter.CRE))
}

GRN_switch_search <- function(GRN.data, seurat.object, cores = 6, enhancers, CRE.score.cutoff=0.35, genes.use = NULL, Pe.to.G.cutoff = 0.05, Pe.to.G.p.value.cutoff = 0.05){
  require(Seurat)
  require(cicero)
  require(tidyverse)
  require(presto)
  require(Signac)
  require(plyr)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(parallel)
  require(pbmcapply)

  # Entire search space in terms of genes
  if (is.null(genes.use)){
    G.a <- unique(GRN.data$markers_rna$RNA.feature)
  } else {
    G.a <- genes.use
  }

  # Prefilter peaks to G to links
  Pe.to.G <- GRN.data$peaks.to.genes[GRN.data$peaks.to.genes$score >= Pe.to.G.cutoff & GRN.data$peaks.to.genes$pvalue <= Pe.to.G.p.value.cutoff,]
  Pe.to.G.posneg <- GRN.data$peaks.to.genes[(GRN.data$peaks.to.genes$score >= Pe.to.G.cutoff | GRN.data$peaks.to.genes$score <= Pe.to.G.cutoff) & GRN.data$peaks.to.genes$pvalue <= 0.05,]
  RNA.markers.ident <- tibble(GRN.data$markers_rna) %>% filter(RNA.group == GRN.data$ident.of.interest, RNA.auc > 0.5, RNA.padj < 0.05) %>% pull(RNA.feature)

  # Take peak coordinates from Pe to Genes list (they are in metadata columns) and form GRange object for temporary use
  Pe.to.G.Pe.coords <- StringToGRanges(GRN.data$peaks.to.genes$peak)

  # Lapply over the entire gene search space
 output <- pbmclapply(G.a, function(G){
 # output <- lapply(G.a, function(G){
   tryCatch({
      #print(G)
      G.output<-list()
      # if given gene id is not Ensembl gene id convert it to such, code also expects mapping to Ensembl gene id in seurat RNA assay feature metadata
      if (!grepl(pattern = "^ENS.*", x = G) & ("ENSEMBL" %in% names(seurat.object[["RNA"]][[]]))){
        G.ensg <- seurat.object[["RNA"]][[]][G,1]
        
        G.output$ensg <- G.ensg
        # Downstream check 1:
        # Set the gene type either G.tf or G.non.tf
        G.output$type <- ifelse(is.G.tf(G, GRN.data$markers_motifs$gene), "G.tf", "G.non.tf")
  
        # Downstream check 2:
        downstream_AUC_tmp <- filter(GRN.data$TF_RNA_TF_z_pairs, gene == G)
        
        # If G.tf == TRUE, and there is AUC data then let's look for downstream target activity via AUC statistics
        if (G.output$type == "G.tf" & !nrow(downstream_AUC_tmp)==0){
          G.output$downstream_filtered <- G.tf.downstream.targets(G, GRN.data, seurat.object)
          G.output$downstream_AUC <- downstream_AUC_tmp
        } else {
          G.output$downstream_filtered <- NA
          G.output$downstream_AUC <- NA
        }
      
        # Downstream check 3:
        # If G.tf == FALSE, then we can check whether gene is expressed marker in the cells
        G.output$RNA.marker <- is.RNA.marker(G, GRN.data,RNA.padj.cutoff =  0.05, RNA.logFC = .5)
        
        # Upstream check 1:
        # Looking for upstream of G what peaks are associated with G.rna
        # TODO: following line is repeated in the fetch.associated.CREs(), should be rewritten later on
        G.output$Pe.to.G <- is.Pe.to.Grna(G, Pe.to.G.posneg, seurat.object)
  
        # Upstream check 2:
        # CRE processing
        fetched.CRE <- fetch.associated.CREs(GRN.data, CRE.score.cutoff, enhancers, G.output)
        G.output$G.non.promoter.CRE <- fetched.CRE$G.non.promoter.CRE
        G.output$G.promoter.CRE <- fetched.CRE$G.promoter.CRE
      } else {
        G.output$G <- G
      }
    return(G.output)
   },
   warning=function(cond){
     return(cond)
   },
   error=function(cond){
     return(cond)
   }
   )
      #}
}, mc.cores = cores
  )

  names(output) <- G.a

  # Then let's calculate global Footprinting as we should see the G.tf in question having some footprint signal if it is really actively regulating target gene expression
  #all.motifs.with.underscore <- gsub(x = GRN.data$markers_motifs$motif.feature, pattern = "-", replacement = "_")
  # Following does not work well as Footprinting signal cannot be reliably calculated from only small subset of regions
  #tmp <- Footprint(object = seurat.object, genome = BSgenome.Mmusculus.UCSC.mm10, regions = G.downstream.filt2, key="Dlx2foot")

  # Following works better as there is more data to resolve
  #seurat.object <- Footprint(object = seurat.object, genome = BSgenome.Mmusculus.UCSC.mm10, motif.name = all.motifs.with.underscore, in.peaks = TRUE)

  #p2 <- PlotFootprint(tmp, features = c("Dlx2foot"), idents = c("P3 GABAergic neurons", "P1 GABA"))
  #p2 <- PlotFootprint(tmp, features = c("Dlx2_MA0885.1"), idents = c("P3 GABAergic neurons", "P1 GABA"))
  # Here problem with the footprinting is interpretation of the results, and that it is still aggregate over all siten given in the regions argument
  return(output)
}

format_GRN_search_data<-function(seurat.object, ident.of.interest, group.by, cell.pseudotime = NULL, padj.cutoff = 1e-2, CRE.score.cutoff = .25, txdb = TxDb.Mmusculus.UCSC.mm10.ensGene, genes.use = NULL){
  # Assuming now that seurat.object contains assays peaks, RNA, chromvar, also links and motif.match for footprinting

  require(Seurat)
  require(cicero)
  require(tidyverse)
  require(presto)
  require(Signac)
  require(TxDb.Mmusculus.UCSC.mm10.ensGene)
  require(org.Mm.eg.db)
  require(clusterProfiler)

  # Initialize output object
  output <- list()

  # Write here smoothing with pseudo time if pseudotime has been given

  # Presto AUC ----

  # Using presto based WX test to find TFs which are both expressed and which motif are accessible
  print("Running presto::wilcoxauc for RNA assay")
  markers_rna <- presto:::wilcoxauc.Seurat(X = seurat.object, group_by = group.by, assay = 'data', seurat_assay = 'RNA') %>% filter(group==ident.of.interest)
  print("Running presto::wilcoxauc for chromvar assay")
  markers_motifs <- presto:::wilcoxauc.Seurat(X = seurat.object, group_by = group.by, assay = 'data', seurat_assay = 'chromvar') %>% filter(group==ident.of.interest)

  # Formatting resulting objects
  motif.names <- markers_motifs$feature
  colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
  colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
  markers_rna$gene <- markers_rna$RNA.feature
  # ConverMotifID does not work with 10x motifs due to the lacking mapping between motif id and name (both catenated to string)
  #markers_motifs$gene <- ConvertMotifID(seurat.object, id = motif.names)
  # Custom code to extract gene names, remove (var) infor, and change to capital first letter format
  markers_motifs$gene <-  str_to_title(sapply(strsplit(motif.names,"-"), function(x){x[1]})) %>% gsub(pattern = "\\(.*\\)", replacement = "")

  # TFrna to TF z-score ----
  # Call for separate function (from Satija lab pages) for detection of interesting TF and TF-motif pairs
  print("Find TF_RNA and TF_z-scorea associations")
  TF_RNA_TF_z_pairs <- topTFs(celltype = ident.of.interest, padj.cutoff = padj.cutoff, markers_motifs = markers_motifs, markers_rna = markers_rna)

  output$markers_rna <- markers_rna
  output$markers_motifs <- markers_motifs
  output$TF_RNA_TF_z_pairs <- tibble(TF_RNA_TF_z_pairs)
  output$ident.of.interest <- ident.of.interest

  # Cicero ----

  # Run cicero specific for the subset
  DefaultAssay(seurat.object) <- 'peaks'
  # convert to CellDataSet format and make the cicero object
  s.data.cds <- as.cell_data_set(x = seurat.object)
  cds.ident.boolean <- pData(s.data.cds)$ident == ident.of.interest
  s.data.cicero <- make_cicero_cds(s.data.cds[,cds.ident.boolean], reduced_coordinates = reducedDims(s.data.cds)$UMAP[cds.ident.boolean,], k=min(50,ceiling(sum(cds.ident.boolean)*.1)))

  # get the chromosome sizes from the Seurat object
  genome <- seqlengths(seurat.object)

  # convert chromosome sizes to a dataframe
  genome.df <- data.frame("chr" = names(genome), "length" = genome)

  conns <- run_cicero(s.data.cicero, genomic_coords = genome.df, sample_num = 100)
  ccans <- generate_ccans(conns)

  links <- ConnectionsToLinks(conns = conns, ccans = ccans)
  output$cicero_links <- links

  # Analyze through Cicero CRE links
  CRE.data <- output$cicero_links
  filt.CRE.data <- CRE.data[CRE.data$score >= CRE.score.cutoff]
  # Let's annotate CRE data, it is granges object thus let's annotate it as such

  txdb_str_chrom <- keepStandardChromosomes(txdb)
  mm10.promoters<- promoters(txdb_str_chrom, upstream = 2000, downstream = 400, columns = "GENEID")

  # Finding which CREs are in promoter area and which not
  promoter.CRE <- filt.CRE.data[filt.CRE.data %over% mm10.promoters]
  non.promoter.CRE <- filt.CRE.data[filt.CRE.data %outside% mm10.promoters]
  output$promoter.CRE <- promoter.CRE
  output$non.promoter.CRE <- non.promoter.CRE

  # Adding Ensembl gene ids to promoter.CRE Grange
  #mcols(promoter.CRE, level="within")$GENEID <- unlist(mm10.promoters[mm10.promoters %over% promoter.CRE]$GENEID)

  # Calculate peak to gene expression associations ----
  peaks.to.genes <- LinkPeaks(subset(seurat.object, idents = ident.of.interest), peak.assay = "peaks", expression.assay = "RNA", n_sample=200, genes.use=genes.use)
  output$peaks.to.genes <- Links(peaks.to.genes)
  
  output$motifs <- Motifs(seurat.object[["peaks"]])

  # Storing  promoter coordinate info for genes
  output$promoter.coordinates <- mm10.promoters

  return(output)
}

topTFs <- function(celltype, padj.cutoff = 1e-2, markers_rna, markers_motifs) {
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
