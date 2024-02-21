# Set of functions to handle TOBIAS data object (Granges list with specific content)

get_BINDetect_snakemake_results_gr <- function(res_path,parallel=F, mc.cores=NULL, HOCOMOCO=11){
  
  #'@param res_path (str): Path to the folder where TOBIAS BINDetect results are stored
  #'
  #'@returns a named list (Large list) consisting of granges for each result sub folder in @param res_path.
  #'         The list can be conveniently accessed, for example, with out_list$gene_TFBSname
  #'
  #'@example get_BINDetect_snakemake_results("/path/to/TOBIAS_framework/outs/TFBS/")
  #'
  # 'Dependencies'
  require(GenomicRanges)
  require(magrittr)
  require(tidyverse)
  require(parallel)
  
  # Reject the .txt, .pdf, etc. files with regex.
  # Apparently all sub folders are of form gene_TFBSname.n where n \in {1,2,3}
  if (HOCOMOCO==11){
    motif.res.folders <- list.files(res_path, pattern = "(.*\\.H[0-9]{2}MO\\.[A-Z]{1})|(\\.[0-9])")
  } else if (HOCOMOCO==12){
    motif.res.folders <- list.files(res_path, pattern = "*H12CORE.*|*H11MO.*")
  }
  # Drop non-folders from the list
  motif.res.folder.i <- sapply(motif.res.folders,function(d){dir.exists(paste(res_path,d,sep=""))})
  motif.res.folders <- motif.res.folders[motif.res.folder.i]
  
  if (!parallel){
    # The actual loop as described in pseudo
    out_list <- lapply(motif.res.folders, function(name) {
      # Access the sub folder's contents.
      # This should be of form res_path/gene_TFBSname.n/beds/
      overview.file.path <- paste0(res_path, name) %>% paste0("/",name,"_overview.txt")
      
      # A little derail, but apparently the most simple way to name each column in the granges is
      # to convert the bed-file into a column-named data frame.
      # The Granges inherits the column names and thus is can be indexed by column names.
      overview.df <- data.frame(read.table(overview.file.path,header = TRUE))
      #selected.cols.i <- TRUE#!grepl(pattern=".*_log2fc",x=colnames(overview.df))
      overview.df.import <- overview.df#[,selected.cols.i]
      colnames(overview.df.import) <- str_replace(string = colnames(overview.df.import), pattern = "_footprints_bound", replacement = "_bound") %>% str_replace(pattern = "_footprints_score", replacement = "_score")
      # TODO: This could be significantly faster if one first catenates all files and reads it as one pass into df and then GR
      GenomicRanges::makeGRangesFromDataFrame(overview.df.import, keep.extra.columns = TRUE,
                                              seqnames.field = "TFBS_chr",
                                              start.field = "TFBS_start",
                                              end.field = "TFBS_end",
                                              strand.field = "TFBS_strand")
    })
    names(out_list) <- motif.res.folders
  } else {
    require(parallel)
    # The actual loop as described in pseudo
    out_list <- mclapply(motif.res.folders, function(name) {
      # Access the sub folder's contents.
      # This should be of form res_path/gene_TFBSname.n/beds/
      overview.file.path <- paste0(res_path, name) %>% paste0("/",name,"_overview.txt")
      
      # A little derail, but apparently the most simple way to name each column in the granges is
      # to convert the bed-file into a column-named data frame.
      # The Granges inherits the column names and thus is can be indexed by column names.
      overview.df <- data.frame(read.table(overview.file.path,header = TRUE))
      selected.cols.i <- !grepl(pattern=".*_log2fc",x=colnames(overview.df))
      overview.df.import <- overview.df[,selected.cols.i]
      colnames(overview.df.import) <- str_replace(string = colnames(overview.df.import), pattern = "_footprints_bound", replacement = "_bound") %>% str_replace(pattern = "_footprints_score", replacement = "_score")
      GenomicRanges::makeGRangesFromDataFrame(overview.df.import, keep.extra.columns = TRUE,
                                              seqnames.field = "TFBS_chr",
                                              start.field = "TFBS_start",
                                              end.field = "TFBS_end",
                                              strand.field = "TFBS_strand")
    }, mc.cores=mc.cores)
    names(out_list) <- motif.res.folders
  }
  return(out_list)
}

combine.footprint.gr <- function(gr.list){
  # Assumes identical motifs per each gr in gr.list
  
  all.motifs <- get.TF.motifs(gr.list[[1]])
  h<-hash()
  for (m in all.motifs){
    for (g in gr.list){
      if(!has.key(m,h)){
        h[[m]]<-as.data.frame(g[[m]])
      } else {
        h[[m]]<-list(h[[m]],as.data.frame(g[[m]]))
      }
      
    }
  }
  return(h)
}


ConstructBed_TobiasGr <- function(gr,group,TF,file=TRUE, gr.only=FALSE){
  require(valr)
  if(any(length(group)>1 | length(TF)>1)){
    errorCondition("Only one group and TF at the time are supported")
  } else {
    TF.ni <- grep(x=names(gr), pattern=str_to_upper(TF))
    beds.out<-lapply(TF.ni, function(v) {
      gr.tmp <- gr[[v]]
      
      group.ni <- which(colnames(elementMetadata(gr.tmp))==paste(group,"_bound",sep=""))
      score.ni <- which(colnames(elementMetadata(gr.tmp))==paste(group,"_score",sep=""))
      
      gr.out <- gr.tmp[(elementMetadata(gr.tmp)[,group.ni]==1)]
      if (!gr.only){
        bed.out <- gr_to_bed(gr.out)
        bed.out$footprint_score <- elementMetadata(gr.out)[,score.ni]
      } else {
        bed.out <- gr.out
      }
      
      if (is.null(file)){
        return(bed.out)
      } else {
        invisible(bed.out)
        write_tsv(bed.out, file=paste(group,"_",names(gr)[v],"_.bed",sep=""), col_names=FALSE)
      }
    })
    
  }
  if(gr.only){
    suppressWarnings(beds.out<-unlist(as(beds.out, "GRangesList")))
    return(beds.out)
  }
}


newPlot<-function(gr,motif,cond1,cond2,dataset){
  gr.i<- which(names(gr)==motif)
  tmp.gr <- as_tibble(gr[[gr.i]]) %>% filter(name %in% c("protein_coding_promoter", "any_promoter"))
  cond.name <- paste(cond1,"_",cond2,"_log2fc",sep="")
  bound.name <- paste(cond1,"_bound",sep="")
  
  if (!cond.name %in% colnames(tmp.gr)){
    cond1.tmp <- cond1
    cond1 <- cond2
    cond2 <- cond1.tmp
    cond.name <- paste(cond1,"_",cond2,"_log2fc",sep="")
    rev.comp <- TRUE
  }
  
  fp.log2fc <- tmp.gr[,cond.name]
  bound <- tmp.gr[,bound.name]
  
  # Calculate expression data fold change
  DefaultAssay(dataset) <- "RNA"
  RNA.fold <- FoldChange(object = dataset, ident.1 = cond1, ident.2 = cond2,assay="RNA") %>% rownames_to_column(var="gene_id") %>% as_tibble() %>% filter(pct.1 > 0.25 | pct.2 > 0.25) %>% filter(!is.na(avg_log2FC)) 
  
  tmp.gr<-inner_join(tmp.gr,RNA.fold)
  
  tb.to.plot <- select(tmp.gr,all_of(cond.name),all_of(bound.name),"avg_log2FC","gene_name") %>% filter(!if_any(everything(), is.na)) %>% arrange(desc(get({{cond.name}})))
  
  matrix.to.plot <- as.matrix(tb.to.plot[,1:3])
  
  col_fun.1 = colorRamp2(c(min(matrix.to.plot[,1]), mean(matrix.to.plot[,1]), max(matrix.to.plot[,1])), c("blue", "white", "red"))
  h1 <- Heatmap(matrix.to.plot[,1], cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, col = col_fun.1)
  
  h2 <- Heatmap(matrix.to.plot[,2], cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, col=c("white","darkgreen"))
  
  col_fun.3 = colorRamp2(c(min(matrix.to.plot[,3]), mean(matrix.to.plot[,3]), max(matrix.to.plot[,3])), c("blue", "white", "red"))
  h3 <- Heatmap(matrix.to.plot[,3], cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, col=col_fun.3)

  h1 + h3 + h2
}

get.conditions <- function(tobias.gr){
  all.metadata <- lapply(tobias.gr,function(f){colnames(elementMetadata(f))})
  all.equal <- all(sapply(all.metadata, identical, all.metadata[[1]]))
  if (all.equal){
    metadata.colnames <- colnames(elementMetadata(tobias.gr[[1]]))
    conditions <- str_extract(string=metadata.colnames, pattern=".*_bound") %>% str_remove(pattern = "_bound") %>% na.omit()
  } else {
    errorCondition("Non-identical metadata columns in TOBIAS gr")
  }
return(as.character(conditions))
}

get.TF.motifs <- function(tobias.gr, remove_repeat=FALSE){
  if (remove_repeat){
    motif.names <- str_remove(string=names(tobias.gr), pattern = "_.*")
    return(motif.names)
  } else {
    return(names(tobias.gr))
  }
}

get.footprints <- function(tobias.gr, conditions, TF.motif, gr.filter=NULL, binary=TRUE){
  require(Signac)
  tob.sub <- tobias.gr[[TF.motif]]
  if (!is.null(gr.filter)){
    if (!class(gr.filter)=="GRanges"){gr.filter <- StringToGRanges(gr.filter)}
    tob.sub <- tob.sub[tob.sub %over% gr.filter]
  }
  if (binary & length(tob.sub)>0){
    bound.info <- as.data.frame(elementMetadata(tob.sub)[,paste(conditions,"_bound",sep="")])
    out <- cbind(coordinate=GRangesToString(tob.sub), bound.info)
    colnames(out) <- c("coordinate",paste(conditions,"_bound",sep=""))
    return(out)
  } else if (!binary & length(tob.sub)>0) {
    bound.info <- as.data.frame(elementMetadata(tob.sub)[,paste(conditions,"_score",sep="")])
    out <- cbind(coordinate=GRangesToString(tob.sub), bound.info)
    colnames(out) <- c("coordinate",paste(conditions,"_score",sep=""))
    return(out)
  } else {
    out <- as.data.frame(matrix(nrow=1, ncol=length(conditions), data = NA))
    if (binary){
      colnames(out) <- paste(conditions,"_bound",sep="")
    } else {
      colnames(out) <- paste(conditions,"_score",sep="")
    }
    return(out)
  }
}

get.fp.logfc <- function(tobias.gr, TF.motif, cond1, cond2, gr.filter){
  require(plyranges)
  tmp.gr <- tobias.gr[[TF.motif]] %>% plyranges::filter_by_overlaps(StringToGRanges(gr.filter))
  tmp.gr.metadata.names <- colnames(elementMetadata(tmp.gr))
  tmp.gr.metadata<-elementMetadata(tmp.gr)
  rev<-FALSE
  if (any(tmp.gr.metadata.names==paste(cond1,cond2,"log2fc",sep="_"))){
    col.ni <- which(tmp.gr.metadata.names==paste(cond1,cond2,"log2fc",sep="_"))
  } else {
    col.ni <- which(tmp.gr.metadata.names==paste(cond2,cond1,"log2fc",sep="_"))
    rev <- TRUE
  }
  
  tmp.out <- unlist(tmp.gr.metadata[col.ni])
  
  if (!rev){
    tmp.out<-tmp.out*-1
  }
  
  if (length(tmp.out)>1){
    return(log2(mean(2^tmp.out)))
  } else {
    return(tmp.out)
    }
}

get.gene.exp.log2fc <- function(dataset, motif,cond1, cond2, group.by, H.metadata,min.exp=0.25){
  gene.id <- H.metadata$ensg_id[match(motif,H.metadata$name)]
  TF.expression.data.avg <- exp(AverageExpression(dataset, assays = "RNA", features = gene.id, group.by = group.by)[[1]][,c(cond1,cond2)])
  tmp.gene.stats <- quantile(dataset@assays$RNA@data[gene.id,],min.exp)
  if (any(TF.expression.data.avg==0) | all(TF.expression.data.avg<tmp.gene.stats)){
    return(NA)
  } else {
    return(log2(TF.expression.data.avg[1]/TF.expression.data.avg[2]))
  }
}


get.gene.exp <- function(dataset, motif,conditions,group.by, H.metadata){
  gene.id <- H.metadata$ensg_id[match(motif,H.metadata$name)]
  TF.expression.data.avg <- AverageExpression(dataset, assays = "RNA", features = gene.id, group.by = group.by)[[1]][,conditions]
  #tmp.gene.stats <- quantile(dataset@assays$RNA@data[gene.id,],min.exp)
  #  all(TF.expression.data.avg<tmp.gene.stats)){
  if (all(TF.expression.data.avg==0)){ 
    return(NA)
  } else {
    return(mean(TF.expression.data.avg,na.rm=T))
  }
}

condence.footprints <- function(footprints.df, binary){
  if (binary){footprint.col.i <- str_detect(string = colnames(footprints.df), pattern = ".*_bound")} else {footprint.col.i <- str_detect(string = colnames(footprints.df), pattern = ".*_score")}
  if (binary & !any(is.na(footprints.df))){
    if (is.null(dim(footprints.df[,footprint.col.i]))){
      tmp.1 <- as.data.frame(sum(footprints.df[,footprint.col.i]))
      colnames(tmp.1) <- colnames(footprints.df)[footprint.col.i]
      tmp.1[tmp.1>0] <-1
    } else {
      tmp.1 <- apply(footprints.df[,footprint.col.i],2,sum)
      tmp.1[tmp.1>0] <-1
    }
    return(tmp.1)
  } else if (!binary & !any(is.na(footprints.df))) {
    if (is.null(dim(footprints.df[,footprint.col.i]))){
      tmp.1 <- mean(footprints.df[,footprint.col.i])
      colnames(tmp.1) <- colnames(footprints.df)[footprint.col.i]
    } else {
      tmp.1 <- apply(footprints.df[,footprint.col.i],2,max)
    }
    return(tmp.1)
  } else {
    return(footprints.df[,footprint.col.i])
  }
}


formFootprintMatrix.overConditions <- function(tobias.gr, conditions, gr.filter, binary=TRUE, na.omit=FALSE){
  all.motifs <- get.TF.motifs(tobias.gr)
  out<-as.data.frame(sapply(all.motifs, function(m){
      tmp.footprints <- get.footprints(tobias.gr, conditions = conditions, TF.motif=m, gr.filter=gr.filter,binary=binary)
      unlist(condence.footprints(tmp.footprints,binary=binary))
    })) 
  rownames(out) <- conditions
  if (na.omit){out<- out[,apply(out,2,function(t){!all(is.na(t))})]}
  return(list(footprint.matrix=as.data.frame(out),feature=gr.filter))
}

readHOCOMOCO.metadata <- function(path){
  require(jsonlite)
  require(dplyr)
  lines <- readLines(path)
  lines <- lapply(lines, fromJSON)
  lines <- lapply(lines, unlist)
  x <- bind_rows(lines)
  return(x)
}

# getTF.gene_symbol <- function(TF.motif, HOCOMOCO.metadata, H.version="H12"){
#   if (H.version=="H12"){gene_symbol <- dplyr::filter(HOCOMOCO.metadata, name %in% TF.motif) %>% pull(masterlist_info.species.MOUSE.gene_symbol)} else {gene_symbol <- dplyr::filter(HOCOMOCO.metadata, Model %in% TF.motif) %>% pull(`Transcription factor`)}
#   return(gene_symbol)
# }

getTF.gene_symbol <- function(TF.motif, HOCOMOCO.metadata, H.version="H12"){
  if (H.version=="H12"){
    #gene_symbol <- dplyr::filter(HOCOMOCO.metadata, name %in% TF.motif) %>% pull(ensg_id)
    tmp.ni <- match(TF.motif, HOCOMOCO.metadata$name)
    gene_symbol <- HOCOMOCO.metadata$masterlist_info.species.MOUSE.gene_symbol[tmp.ni]
  } #else {gene_symbol <- dplyr::filter(HOCOMOCO.metadata, Model %in% TF.motif) %>% pull(`Transcription factor`)}
  return(gene_symbol)
}

getTF.id <- function(TF.motif, HOCOMOCO.metadata, H.version="H12"){
  if (H.version=="H12"){
    #gene_symbol <- dplyr::filter(HOCOMOCO.metadata, name %in% TF.motif) %>% pull(ensg_id)
    tmp.ni <- match(TF.motif, HOCOMOCO.metadata$name)
    gene_id <- HOCOMOCO.metadata$ensg_id[tmp.ni]
  } #else {gene_symbol <- dplyr::filter(HOCOMOCO.metadata, Model %in% TF.motif) %>% pull(`Transcription factor`)}
  return(gene_id)
}



plotTFfootprint.heatmap <- function(footprint.matrix, filter.unbound=FALSE, with.expression=FALSE, filter.no.expressed.TF=FALSE, Seurat.dataset=NULL, HOCOMOCO.metadata=NULL, H.version="H12", group.by="rv2.lineage"){
  require(ComplexHeatmap)
  require(circlize)
  require(patchwork)
  require(factoextra)
  
  feature <- footprint.matrix$feature
  footprint.matrix<-footprint.matrix$footprint.matrix
  
  conditions <- str_remove(string = rownames(footprint.matrix), pattern = "_bound") %>% str_remove(pattern = "_score")
  rownames(footprint.matrix) <- conditions
  if (H.version=="H12"){colnames(footprint.matrix) <- str_remove(string=colnames(footprint.matrix), pattern = "_.*")} else {colnames(footprint.matrix) <- str_remove(string = colnames(footprint.matrix), pattern = "^.*\\.[:alpha:]+_")}
  if (filter.unbound){keep.i <- which(colSums(footprint.matrix)>0); footprint.matrix<-footprint.matrix[,keep.i]}

  if (with.expression & !is.null(Seurat.dataset) & !is.null(HOCOMOCO.metadata)){
    TF.motif.names <- colnames(footprint.matrix)
    TF.ensg <- getTF.id(TF.motif = TF.motif.names, HOCOMOCO.metadata = HOCOMOCO.metadata, H.version = H.version)
    
    DefaultAssay(Seurat.dataset) <- "RNA"
    TF.expression.data.avg <- AverageExpression(Seurat.dataset, assays = "RNA", features = na.omit(TF.ensg), group.by = group.by)[[1]][,conditions]
    # Repeat some TF expression to match motifs repeats
    TF.expression.data.avg<-TF.expression.data.avg[TF.ensg,]
    
    # Switch expression data to have motif names
    rownames(TF.expression.data.avg) <- TF.motif.names

    # Filter non-expressed TFs
    if (filter.no.expressed.TF){
      TF.i.2 <- rowMeans(TF.expression.data.avg)>quantile(TF.expression.data.avg, 0.5)
      footprint.matrix <- footprint.matrix[,TF.i.2]
      TF.expression.data.avg<-TF.expression.data.avg[TF.i.2,]
    }
    # Drop genes having NA expression
    TF.expression.data.avg.scaled<-t(scale(t(TF.expression.data.avg)))
    #color.func.exp <- colorRamp2(c(min(TF.expression.data.avg.scaled), mean(apply(TF.expression.data.avg.scaled,1,mean)), max(TF.expression.data.avg.scaled)), c("white", "white", "red"))
  }
  
  #feat.accessibility.data.avg <- AverageExpression(Seurat.dataset, assays = "peaks", features = feature, group.by = group.by)[[1]][,conditions]
  
  footprint.matrix.scaled <- scale(t(footprint.matrix))
  
  combined.matrix <- cbind(TF.expression.data.avg.scaled,footprint.matrix.scaled)

  na.ni <- which(rowSums(is.na(combined.matrix))>0)
  if (length(na.ni)>0){
    combined.matrix<- combined.matrix[-na.ni,]
  }
  fviz_nbclust.res <- fviz_nbclust(combined.matrix, FUN=kmeans, method="silhouette")
  opt.kmeans <-  which.max(fviz_nbclust.res$data$y)
  
  TF.split <- kmeans(combined.matrix,centers=opt.kmeans, iter.max = 100)
  color.func.foot <- colorRamp2(c(min(combined.matrix), mean(apply(combined.matrix,2,mean)), max(combined.matrix)), c("blue", "white", "red"))
    
  p.foot.exp <- Heatmap(combined.matrix, cluster_columns = FALSE, cluster_rows = TRUE, row_names_gp = gpar(fontsize = 6), col = color.func.foot, split=TF.split$cluster, column_split = factor(c(rep("expression",5),rep("fp",5)), levels=c("expression","fp")), column_gap = unit(5, "mm"))
  
  return(p.foot.exp)
    
  # column_ha = HeatmapAnnotation(avg.acc = anno_barplot(feat.accessibility.data.avg), col=list(avg.acc="darkgreen"))
  # 
  # #foot.dist <- dist(t(footprint.matrix), method="euclidean")
  # #foot.hclust <- hclust(foot.dist, method = "complete")
  # 
  # if (!any(footprint.matrix>0 & footprint.matrix<1)){color.func.foot <- c("white","darkgreen")} else {color.func.foot <- colorRamp2(c(min(footprint.matrix.scaled), mean(apply(footprint.matrix.scaled,2,mean)), max(footprint.matrix.scaled)), c("white", "yellow", "darkgreen"))}
  # p.foot <- Heatmap(footprint.matrix.scaled, cluster_columns = FALSE, cluster_rows = TRUE, row_names_gp = gpar(fontsize = 6), col = color.func.foot, split=TF.split$cluster, top_annotation = column_ha)
  # 
  # if (with.expression){
  #   p.exp <- Heatmap(TF.expression.data.avg.scaled, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize = 6), col = color.func.exp,split=TF.split$cluster)
  #   p.final <- p.foot + p.exp + plot_layout(ncol=2)
  #   return(p.final)
  # } else {
  #   return(p.foot)
  # }
    
  
}


plotFootprintDotplot.rV2.comp <- function(tobias.fp.gr, gr.filter, TF.motifs, dataset, H.metadata=H12.metadata, fp.pos.thr=1, fp.neg.thr=-1, exp.pos.thr=0.5, exp.neg.thr=-0.5, min.exp=.1, verbose=T, mc.cores=9){
  # Arguments:
  # tobias.fp.gr = TOBIAS dataobject a such from qs file
  # gr.filter = string of genomic coordinates of interest chr-start-end
  # TF.motifs = tibble with motif names and expression.class
  # dataset = rV2 Seurat object with groups merged into rv2.lineage_re metadata variable
  # H.metadata = Hocomoco metadata information from H12_metadata_mod.qs
  # Could be optimized a lot by building everything in one loop over TFs, but makes no sense to do that now
  
  require(parallel)
  DefaultAssay(dataset) <- "RNA"
  
  # Preparing data
  TF.motifs.rep <- paste(TF.motifs$TF.motif,TF.motifs$TF.motif,sep="_")
  
  if (verbose){print("Preparing data from PRO/CO transition")}
  
  tmp.CO <- mclapply(TF.motifs.rep,function(tf){
    get.fp.logfc(tobias.gr = rV2.groups.tobias.h12.gr.dr, TF.motif = tf, cond1="PRO1_2",cond2="CO1_2", gr.filter = gr.filter)
    }, mc.cores=mc.cores)
  names(tmp.CO) <- TF.motifs.rep
  tmp.CO <- unlist(tmp.CO)

  CO.tb <- tibble(group="PRO/CO",motif=names(tmp.CO),fp.log2fc=tmp.CO)
  CO.tb$motif <- str_remove(CO.tb$motif, pattern = "_.*")
  CO.tb$TF.exp.log2fc <- unlist(mclapply(CO.tb$motif,function(m){
    get.gene.exp.log2fc(dataset=dataset, motif=m,cond1="PRO1_2",cond2="CO1_2",group.by = "rv2.lineage_re", H.metadata = H.metadata)
  },mc.cores=mc.cores))

  
  if (verbose){print("Preparing data from CO/GA transition")}
  
  tmp.GA <- mclapply(TF.motifs.rep,function(tf){
    get.fp.logfc(tobias.gr = rV2.groups.tobias.h12.gr.dr, TF.motif = tf, cond1="CO1_2",cond2="GA1_2", gr.filter = gr.filter)
  }, mc.cores=mc.cores)
  names(tmp.GA) <- TF.motifs.rep
  tmp.GA <- unlist(tmp.GA)

  GA.tb <- tibble(group="CO/GA",motif=names(tmp.GA),fp.log2fc=tmp.GA)
  GA.tb$motif <- str_remove(GA.tb$motif, pattern = "_.*")
  GA.tb$TF.exp.log2fc <- unlist(mclapply(GA.tb$motif,function(m){
    get.gene.exp.log2fc(dataset=dataset, motif=m,cond1="CO1_2",cond2="GA1_2",group.by = "rv2.lineage_re", H.metadata = H.metadata)
  },mc.cores=mc.cores))

  
  if (verbose){print("Preparing data from CO/GL transition")}
  
  tmp.GL <- mclapply(TF.motifs.rep,function(tf){
    get.fp.logfc(tobias.gr = rV2.groups.tobias.h12.gr.dr, TF.motif = tf, cond1="CO1_2",cond2="GL1_2", gr.filter = gr.filter)
  }, mc.cores=mc.cores)
  names(tmp.GL) <- TF.motifs.rep
  tmp.GL <- unlist(tmp.GL)

  GL.tb <- tibble(group="CO/GL",motif=names(tmp.GL),fp.log2fc=tmp.GL)
  GL.tb$motif <- str_remove(GL.tb$motif, pattern = "_.*")
  
  GL.tb$TF.exp.log2fc <- unlist(mclapply(GL.tb$motif,function(m){
    get.gene.exp.log2fc(dataset=dataset, motif=m,cond1="CO1_2",cond2="GL1_2",group.by = "rv2.lineage_re", H.metadata = H.metadata)
  },mc.cores=mc.cores))

  GA.GL.tb <- full_join(GA.tb,GL.tb)
  GA.GL.CO.tb <- full_join(GA.GL.tb,CO.tb)
  
  # GA.GL.CO.tb$fp.log2fc <- GA.GL.CO.tb$fp.log2fc*-1
  
  GA.GL.CO.tb<- left_join(GA.GL.CO.tb, TF.motifs, by=c("motif"="TF.motif"))
  GA.GL.CO.tb$TF.exp.log2fc <- GA.GL.CO.tb$TF.exp.log2fc*-1
  
  if (verbose){print("Generating plot")}
  
  # Plotting
  data.2.plot <- filter(GA.GL.CO.tb, (fp.log2fc > fp.pos.thr | fp.log2fc < fp.neg.thr) & (TF.exp.log2fc > exp.pos.thr | TF.exp.log2fc < exp.neg.thr) & !(expression.class=="Uncorrelated"))
  data.2.plot$direction <- factor(x = ifelse(data.2.plot$fp.log2fc>0,"Positive","Negative"), levels=c("Positive","Negative"))
  data.2.plot$group <- factor(data.2.plot$group, levels=c("PRO/CO","CO/GA","CO/GL"))
  data.2.plot$motif.with.class <- ifelse(data.2.plot$expression.class=="Correlated",paste(data.2.plot$motif,"(C)",sep=" "), paste(data.2.plot$motif,"(S)",sep=" "))

    p.1 <- ggplot() + 
    geom_point(aes(x = group, y=motif.with.class, size=abs(fp.log2fc), fill=TF.exp.log2fc, shape=direction,group=direction),filter(data.2.plot,group=="PRO/CO")) + scale_shape_manual(values=c(24,25)) +
    geom_point(aes(x = group, y=motif.with.class, size=abs(fp.log2fc), fill=TF.exp.log2fc, shape=direction,group=direction),filter(data.2.plot,group=="CO/GA")) + scale_shape_manual(values=c(24,25)) +
    geom_point(aes(x = group, y=motif.with.class, size=abs(fp.log2fc), fill=TF.exp.log2fc, shape=direction,group=direction),filter(data.2.plot,group=="CO/GL")) + scale_shape_manual(values=c(24,25)) +
    theme_minimal() + scale_fill_gradient2(low="blue", mid="gray",high="red", midpoint=0) + scale_x_discrete(drop = FALSE) + scale_shape_manual(name = "FP change",
                                                                                                                                               labels = c("Positive", "Negative"),
                                                                                                                                               values = c(24,25))
  
  return(p.1)
}


fetchDotPlotData <- function(tobias.gr, TF.motifs, conditions, gr.filter, cons.gr, dataset, mc.cores){
  require(parallel)
  tmp.1 <- mclapply(TF.motifs,function(tf){
    tmp.scores <- get.footprints(tobias.gr = tobias.gr, TF.motif = tf, conditions=conditions, gr.filter = gr.filter, binary = FALSE)
    if (ncol(tmp.scores)==2){
      if (!is.null(cons.gr)){
        tmp.mean.cons <- mean(filter_by_overlaps(cons.gr,StringToGRanges(tmp.scores[,1]))$score, na.rm=T)
        return(c("fp.mean"=mean(tmp.scores[,2]), "cons.mean"=tmp.mean.cons))
      } else {
        return(mean(tmp.scores[,2]))
      }
    } else {
      return(NA)
    }
  }, mc.cores=mc.cores)
  names(tmp.1) <- TF.motifs
  
  if (is.null(cons.gr)){
    tmp.1 <- unlist(tmp.1)
    tmp.2 <- tibble(group=conditions,motif=names(tmp.1),fp=tmp.1)
  } else {
    tmp.2 <- as.data.frame(t(as_tibble(tmp.1))) %>% rownames_to_column("motif")
    colnames(tmp.2) <- c("motif","fp","cons")
    tmp.2$group <- conditions
  }
  
  tmp.2$motif <- str_remove(tmp.2$motif, pattern = "_.*")
  
  tmp.2$TF.exp <- unlist(mclapply(tmp.2$motif,function(m){
    get.gene.exp(dataset=dataset, motif=m,conditions=conditions,group.by = "rv2.lineage_re", H.metadata = H.metadata)
  }, mc.cores = mc.cores))
  
  return(tmp.2)
}


fetch.TF.evidence <- function(tobias.gr, TF.motifs, conditions, gr.filter, cons.gr, dataset, mc.cores){
  require(parallel)
  tmp.1 <- mclapply(TF.motifs,function(tf){
    tmp.scores <- get.footprints(tobias.gr = tobias.gr, TF.motif = tf, conditions=conditions, gr.filter = gr.filter, binary = FALSE)
    if (ncol(tmp.scores)==2){
      if (!is.null(cons.gr)){
        tmp.mean.cons <- mean(filter_by_overlaps(cons.gr,StringToGRanges(tmp.scores[,1]))$score, na.rm=T)
        return(c("fp.mean"=mean(tmp.scores[,2]), "cons.mean"=tmp.mean.cons))
      } else {
        return(mean(tmp.scores[,2]))
      }
    } else {
      return(NA)
    }
  }, mc.cores=mc.cores)
  names(tmp.1) <- TF.motifs
  
  if (is.null(cons.gr)){
    tmp.1 <- unlist(tmp.1)
    tmp.2 <- tibble(group=conditions,motif=names(tmp.1),fp=tmp.1)
  } else {
    tmp.2 <- as.data.frame(t(as_tibble(tmp.1))) %>% rownames_to_column("motif")
    colnames(tmp.2) <- c("motif","fp","cons")
    tmp.2$group <- conditions
  }
  
  tmp.2$motif <- str_remove(tmp.2$motif, pattern = "_.*")
  
  tmp.2$TF.exp <- unlist(mclapply(tmp.2$motif,function(m){
    get.gene.exp(dataset=dataset, motif=m,conditions=conditions,group.by = "rv2.lineage_re", H.metadata = H.metadata)
  }, mc.cores = mc.cores))
  
  return(tmp.2)
}


plotFootprintDotplot.rV2.noncomp <- function(tobias.fp.gr, gr.filter, TF.motifs, dataset, H.metadata=H12.metadata, verbose=F, mc.cores=6, cons.gr=NULL){
  # Arguments:
  # tobias.fp.gr = TOBIAS dataobject a such from qs file
  # gr.filter = string of genomic coordinates of interest chr-start-end
  # TF.motifs = tibble with motif names and expression.class
  # dataset = rV2 Seurat object with groups merged into rv2.lineage_re metadata variable
  # H.metadata = Hocomoco metadata information from H12_metadata_mod.qs
  # Could be optimized a lot by building everything in one loop over TFs, but makes no sense to do that now
  
  DefaultAssay(dataset) <- "RNA"
  
  if (!is.null(cons.gr)){
    cons.subset <- filter_by_overlaps(cons.gr,StringToGRanges(gr.filter))
    }
  
  # Preparing data
  TF.motifs.rep <- paste(TF.motifs$TF.motif,TF.motifs$TF.motif,sep="_")
  
  if (verbose){print("Preparing data for PRO")}
  tmp.PRO <- fetchDotPlotData(tobias.gr=rV2.groups.tobias.h12.gr.dr, TF.motifs=TF.motifs.rep, conditions="PRO1_2", gr.filter=gr.filter, cons.gr=cons.subset, dataset=dataset, mc.cores=mc.cores)
  
  if (verbose){print("Preparing data for CO")}
  tmp.CO <- fetchDotPlotData(tobias.gr=rV2.groups.tobias.h12.gr.dr, TF.motifs=TF.motifs.rep, conditions="CO1_2", gr.filter=gr.filter, cons.gr=cons.subset, dataset=dataset, mc.cores=mc.cores)
  
  if (verbose){print("Preparing data for GA")}
  tmp.GA <- fetchDotPlotData(tobias.gr=rV2.groups.tobias.h12.gr.dr, TF.motifs=TF.motifs.rep, conditions="GA1_2", gr.filter=gr.filter, cons.gr=cons.subset, dataset=dataset, mc.cores=mc.cores)
  
  if (verbose){print("Preparing data for GL")}
  tmp.GL <- fetchDotPlotData(tobias.gr=rV2.groups.tobias.h12.gr.dr, TF.motifs=TF.motifs.rep, conditions="GL1_2", gr.filter=gr.filter, cons.gr=cons.subset, dataset=dataset, mc.cores=mc.cores)
  
  dotplot.data <- full_join(tmp.PRO,tmp.CO)
  dotplot.data <- full_join(dotplot.data,tmp.GA)
  dotplot.data <- full_join(dotplot.data,tmp.GL)
  dotplot.data<- left_join(dotplot.data, TF.motifs, by=c("motif"="TF.motif"))
  
  exp.thr <- mean(dotplot.data$TF.exp,na.rm=T)
  fp.thr <- quantile(dotplot.data$fp,na.rm=T, .15)
  dotplot.data$fp[is.na(dotplot.data$fp)] <- 0
  
  TF.exp.filt.li <- dotplot.data %>% group_by(motif) %>% summarise(TF.filt=any(TF.exp>exp.thr))
  TF.fp.filt.li <- dotplot.data %>% group_by(motif) %>% summarise(TF.fp=any(fp>fp.thr))
  
  dotplot.data <- left_join(dotplot.data,TF.exp.filt.li)
  dotplot.data <- left_join(dotplot.data, TF.fp.filt.li)
  
  dotplot.data$fp[dotplot.data$fp<fp.thr] <- NA
  zero.epression.li <- dotplot.data$TF.exp<exp.thr
  dotplot.data$TF.exp<-log10(dotplot.data$TF.exp)
  dotplot.data$TF.exp[zero.epression.li] <- NA
  
  if (verbose){print("Generating plot")}
  
  # Plotting
  data.2.plot <- filter(dotplot.data,!(expression.class=="Uncorrelated") & TF.filt & TF.fp)
  data.2.plot$group <- factor(data.2.plot$group, levels=c("PRO1_2","CO1_2","GA1_2","GL1_2","cons"))
  data.2.plot$motif.with.class <- ifelse(data.2.plot$expression.class=="Correlated",paste(data.2.plot$motif,"(C)",sep=" "), paste(data.2.plot$motif,"(S)",sep=" "))
  
  p.1 <- ggplot() + 
    geom_point(aes(x = group, y=motif.with.class, size=fp, fill=TF.exp),filter(data.2.plot,group=="PRO1_2"), shape=21) +
    geom_point(aes(x = group, y=motif.with.class, size=fp, fill=TF.exp),filter(data.2.plot,group=="CO1_2"), shape=21)  +
    geom_point(aes(x = group, y=motif.with.class, size=fp, fill=TF.exp),filter(data.2.plot,group=="GA1_2"), shape=21) + 
    geom_point(aes(x = group, y=motif.with.class, size=fp, fill=TF.exp),filter(data.2.plot,group=="GL1_2"), shape=21) + 
    scale_fill_gradient(low="white",high="darkred") +
    new_scale_fill() + 
    geom_point(aes(x="cons", y=motif.with.class, fill=cons),filter(data.2.plot,group=="PRO1_2"),shape=22, size=15) +
    scale_fill_gradient2(low="white",mid = "yellow",high="darkgreen", midpoint = .5) +
    theme_minimal()  + scale_x_discrete(drop = FALSE)
  
  return(list(plot=p.1,data2plot=data.2.plot, all.data=dotplot.data))
}

writeBED.SQL.data <- function(df, filename, meta_column){
  tmp <- select(df, seqnames,start,end,eval(meta_column))
  if (is.numeric(unlist(as.vector(tmp[,4])))){
    write_tsv(tmp,file=paste(filename,".bedgraph",sep=""), col_names = FALSE)
  } else {
    write_tsv(tmp,file=paste(filename,".bed",sep=""), col_names = FALSE)
  }
}

plotHorizDotplot <- function(dbname = "~/Workspace/TOBIAS.dr.h12.sqlite", feature.coords, exp.thr=1.5, mean_cons_thr=.5){
  con.obj <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbname)
  tobias.table <- tbl(con.obj, "tobias")
  exp.table <- tbl(con.obj, "exp_scaled")
  acc.table <- tbl(con.obj, "acc")
  table.tmp.1 <- dplyr::filter(tobias.table, features==feature.coords) %>% left_join(exp.table) %>% left_join(acc.table, by=c("features"="features"))
  
  table.tmp.2 <- table.tmp.1 %>% collect()
  
  table.tmp.2 <- table.tmp.2 %>% filter((abs(PRO1_2.x)>exp.thr | abs(CO1_2.x)>exp.thr | abs(GA1_2.x)>exp.thr) & (PRO1_2_bound==1 | CO1_2_bound==1 | GA1_2_bound==1) & mean_cons>mean_cons_thr) %>% arrange(start)
  
  # Replace cases where TF motif repeats and overlaps itself with max values?.
  
  table.tmp.2$TFBS_name_comb <- paste(str_remove(string = table.tmp.2$TFBS_name, pattern = "_.*")," (",table.tmp.2$start,"-",table.tmp.2$end,")",sep="")
  table.tmp.2$TFBS_name_comb <- factor(table.tmp.2$TFBS_name_comb, levels = table.tmp.2$TFBS_name_comb)
  
  p1 <- ggplot(table.tmp.2) + 
    geom_point(aes(x = TFBS_name_comb, fill=PRO1_2.x, y="4",
                   shape=as.character(PRO1_2_bound), size=PRO1_2_score)) +
    geom_point(aes(x = TFBS_name_comb, fill=CO1_2.x, y="3",
                   shape=as.character(CO1_2_bound), size=CO1_2_score)) + 
    geom_point(aes(x = TFBS_name_comb, fill=GA1_2.x, y="2",
                   shape=as.character(GA1_2_bound), size=GA1_2_score)) +
    geom_point(aes(x = TFBS_name_comb, fill=GL1_2.x, y="1",
                   shape=as.character(GL1_2_bound), size=GL1_2_score)) +
    scale_fill_gradient2(low="blue", mid="gray",high="red", midpoint = 0) + theme_minimal() + theme(axis.text.x = element_text(size = 12, angle = 90)) + ylab("Cell group") + scale_shape_manual(values=c(4,21) , guide = "none") + scale_y_discrete("Cell group", labels=c("4"="PRO1_2","3"="CO1_2","2"="GA1_2","1"="GL1_2")) + labs(fill="Expression z-score", size="Footprint score")

  DBI::dbDisconnect(con.obj)
  return(p1)
  }


plotHorizDotplot_v2 <- function(dbname = "~/Workspace/TOBIAS.dr.h12_2.sqlite", feature.coords, exp.thr=1.2, mean_cons_thr=.5, max.exp, max.acc, max.fp){
  con.obj <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbname)
  tobias.table <- tbl(con.obj, "tobias")
  exp.table <- tbl(con.obj, "exp")
  acc.table <- tbl(con.obj, "acc")
  table.tmp.1 <- dplyr::filter(tobias.table, features==feature.coords) %>% left_join(exp.table) %>% left_join(acc.table, by=c("features"="features"))
  
  table.tmp.2 <- table.tmp.1 %>% collect()
  
  table.tmp.2 <- table.tmp.2 %>% filter((abs(PRO1_2.x)>exp.thr | abs(CO1_2.x)>exp.thr | abs(GA1_2.x)>exp.thr) & (PRO1_2_bound==1 | CO1_2_bound==1 | GA1_2_bound==1) & mean_cons>mean_cons_thr) %>% arrange(start)
  
  # Replace cases where TF motif repeats and overlaps itself with max values?.
  if (nrow(table.tmp.2)==0){return(NA)}
  table.tmp.2$TFBS_name_comb <- paste(str_replace(string = table.tmp.2$TFBS_name, pattern = "([:alpha:])+_.*", "\\1")," (",table.tmp.2$start,"-",table.tmp.2$end,")",sep="")
  table.tmp.2$TFBS_name_comb <- factor(table.tmp.2$TFBS_name_comb, levels = table.tmp.2$TFBS_name_comb)
  
  p1 <- ggplot(table.tmp.2) + 
    geom_point(aes(x = TFBS_name_comb, fill=PRO1_2_score, y="4",
                   shape=as.character(PRO1_2_bound), size=log1p(PRO1_2.x))) +
    geom_point(aes(x = TFBS_name_comb, fill=CO1_2_score, y="3",
                   shape=as.character(CO1_2_bound), size=log1p(CO1_2.x))) + 
    geom_point(aes(x = TFBS_name_comb, fill=GA1_2_score, y="2",
                   shape=as.character(GA1_2_bound), size=log1p(GA1_2.x))) +
    geom_point(aes(x = TFBS_name_comb, fill=GL1_2_score, y="1",
                   shape=as.character(GL1_2_bound), size=log1p(GL1_2.x))) +
    scale_fill_gradient2(low="blue", mid="gray",high="red", midpoint = 0) + theme_minimal() + theme(axis.text.x = element_text(size = 12, angle = 90)) + ylab("Cell group") + scale_shape_manual(values=c(16,21) , guide = "none") + scale_y_discrete("Cell group", labels=c("4"="PRO1_2","3"="CO1_2","2"="GA1_2","1"="GL1_2")) + labs(fill="Footprint score", size="Expression log1p") + xlab("TF-motif") +
    scale_size_continuous(limit=c(0,max.exp))
  
  acc.tmp.dat <- select(table.tmp.2, ends_with(".y")) %>% pivot_longer(everything()) %>% distinct() %>% filter(name %in% c("PRO1_2.y","CO1_2.y","GA1_2.y","GL1_2.y"))
  acc.tmp.dat$name <- acc.tmp.dat$name %>% str_remove(pattern = "\\.y")
  acc.tmp.dat$name <- factor(acc.tmp.dat$name, levels=rev(c("PRO1_2","CO1_2","GA1_2","GL1_2")))
  p2 <- ggplot(acc.tmp.dat, aes(x="Acc",y=name,fill=value)) + geom_point(shape=21, size=5) + theme_minimal() + ylab("Cell group") + xlab("") + scale_fill_viridis(option = "magma", limit=c(0,max.acc)) + labs(fill="Accessibility")
  DBI::dbDisconnect(con.obj)
  return(list(p1=p1,p2=p2))
}

plotHorizDotplot_v3 <- function(dbname = "~/Workspace/TOBIAS.dr.h12_2.sqlite", feature.coords, exp.thr=1.2, mean_cons_thr=.5, igv=NULL,cons.filt=NULL){
  con.obj <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbname)
  tobias.table <- tbl(con.obj, "tobias")
  exp.table <- tbl(con.obj, "exp")
  acc.table <- tbl(con.obj, "acc")
  table.tmp.1 <- dplyr::filter(tobias.table, features==feature.coords) %>% left_join(exp.table) %>% left_join(acc.table, by=c("features"="features"))
  
  table.tmp.2 <- table.tmp.1 %>% collect()
  
  table.tmp.2 <- table.tmp.2 %>% filter((abs(PRO1_2.x)>exp.thr | abs(CO1_2.x)>exp.thr | abs(GA1_2.x)>exp.thr) & (PRO1_2_bound==1 | CO1_2_bound==1 | GA1_2_bound==1) & mean_cons>mean_cons_thr) %>% arrange(start)
  
  # Replace cases where TF motif repeats and overlaps itself with max values?.
  if (nrow(table.tmp.2)==0){return(NA)}
  table.tmp.2$TFBS_name_comb <- paste(str_replace(string = table.tmp.2$TFBS_name, pattern = "([:alpha:])+_.*", "\\1")," (",table.tmp.2$start,"-",table.tmp.2$end,")",sep="")
  table.tmp.2$TFBS_name_comb <- factor(table.tmp.2$TFBS_name_comb, levels = table.tmp.2$TFBS_name_comb)
  
  # Create dotplot
  p1 <- ggplot(table.tmp.2) + 
    geom_point(aes(x = TFBS_name_comb, fill=PRO1_2_score, y="4",
                   shape=as.character(PRO1_2_bound), size=log1p(PRO1_2.x))) +
    geom_point(aes(x = TFBS_name_comb, fill=CO1_2_score, y="3",
                   shape=as.character(CO1_2_bound), size=log1p(CO1_2.x))) + 
    geom_point(aes(x = TFBS_name_comb, fill=GA1_2_score, y="2",
                   shape=as.character(GA1_2_bound), size=log1p(GA1_2.x))) +
    geom_point(aes(x = TFBS_name_comb, fill=GL1_2_score, y="1",
                   shape=as.character(GL1_2_bound), size=log1p(GL1_2.x))) +
    scale_fill_gradient2(low="blue", mid="gray",high="red", midpoint = 0) + theme_minimal() + theme(axis.text.x = element_text(size = 12, angle = 90)) + ylab("Cell group") + scale_shape_manual(values=c(16,21) , guide = "none") + scale_y_discrete("Cell group", labels=c("4"="PRO1_2","3"="CO1_2","2"="GA1_2","1"="GL1_2")) + labs(fill="Footprint score", size="Expression log1p") + xlab("TF-motif")
  
  acc.tmp.dat <- select(table.tmp.2, ends_with(".y")) %>% pivot_longer(everything()) %>% distinct() %>% filter(name %in% c("PRO1_2.y","CO1_2.y","GA1_2.y","GL1_2.y"))
  acc.tmp.dat$name <- acc.tmp.dat$name %>% str_remove(pattern = "\\.y")
  acc.tmp.dat$name <- factor(acc.tmp.dat$name, levels=rev(c("PRO1_2","CO1_2","GA1_2","GL1_2")))
  
  # Create accessibility side plot
  p2 <- ggplot(acc.tmp.dat, aes(x="Acc",y=name,fill=value)) + geom_point(shape=21, size=5) + theme_minimal() + ylab("Cell group") + xlab("") + scale_fill_viridis(option = "magma") + labs(fill="Accessibility")
  
  # Create igvR plot 
  setBrowserWindowTitle(igv, "mm10 TOBIAS")
  setGenome(igv, "mm10")
  showGenomicRegion(igv, str_replace(string = feature.coords, pattern = "-", replacement = ":"))
  
  current.loc <- StringToGRanges(paste(getGenomicRegion(igv)$chrom, getGenomicRegion(igv)$start, getGenomicRegion(igv)$end, sep="-"))
  cons.filt.subset <- cons.filt[cons.filt %over% current.loc]
  cons_filt <- igvR::GRangesQuantitativeTrack("conservation",cons.filt.subset, min=0.5, autoscale = FALSE, max=1)
  displayTrack(igv, cons_filt)
  
  TFBS_to_igvR <- data.frame(table.tmp.2 %>% select(seqnames,start,end))
  TFBS_to_igvR.trck <- DataFrameAnnotationTrack("TFBS", TFBS_to_igvR, color="forestgreen", displayMode="EXPANDED", trackHeight = 200)
  displayTrack(igv, TFBS_to_igvR.trck)
  
  saveToSVG(igv,"igv.tmp.svg")
  
  p3 <- cowplot::ggdraw() + cowplot::draw_image("./igv.tmp.svg")
  
  DBI::dbDisconnect(con.obj)
  return(list(p1=p1,p2=p2,p3=p3))
}

find.maxes <- function(dbname = "~/Workspace/TOBIAS.dr.h12_2.sqlite",features){
  feat.tmp <- pull(features, V2)
  con.obj <- DBI::dbConnect(RSQLite::SQLite(), dbname = "~/Workspace/TOBIAS.dr.h12_2.sqlite")
  tobias.table <- tbl(con.obj, "tobias")
  exp.table <- tbl(con.obj, "exp")
  acc.table <- tbl(con.obj, "acc")
  table.tmp.1 <- dplyr::filter(tobias.table, features %in% feat.tmp) %>% left_join(exp.table) %>% left_join(acc.table, by=c("features"="features")) %>% filter(mean_cons>0.5)
  table.tmp.2 <- table.tmp.1 %>% collect()
  
  max.acc <- max(table.tmp.2 %>% select(ends_with(".y")),na.rm = T)
  max.exp <- log1p(max(table.tmp.2 %>% select(ends_with(".x")), na.rm = T))
  max.fp <- max(table.tmp.2 %>% select(ends_with("_score")), na.rm = T)
  DBI::dbDisconnect(con.obj)
  return(list(max.acc=max.acc, max.exp=max.exp, max.fp=max.fp))
}

