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
    motif.res.folders <- list.files(res_path, pattern = "*H12CORE.*")
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
      selected.cols.i <- !grepl(pattern=".*_log2fc",x=colnames(overview.df))
      overview.df.import <- overview.df[,selected.cols.i]
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
    bound.info <- elementMetadata(tob.sub)[,paste(conditions,"_bound",sep="")]
    out <- cbind(coordinate=GRangesToString(tob.sub), bound.info)
    return(out)
  } else if (!binary & length(tob.sub)>0) {
    bound.info <- elementMetadata(tob.sub)[,paste(conditions,"_score",sep="")]
    out <- cbind(coordinate=GRangesToString(tob.sub), bound.info)
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

condence.footprints <- function(footprints.df, binary){
  if (binary){footprint.col.i <- str_detect(string = colnames(footprints.df), pattern = ".*_bound")} else {footprint.col.i <- str_detect(string = colnames(footprints.df), pattern = ".*_score")}
  if (binary & !any(is.na(footprints.df))){
    tmp.1 <- apply(footprints.df[,footprint.col.i],2,sum)
    tmp.1[tmp.1>0] <-1
    return(tmp.1)
  } else if (!binary & !any(is.na(footprints.df))) {
    tmp.1 <- apply(footprints.df[,footprint.col.i],2,mean)
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
  if (na.omit){out<- out[,apply(out,2,function(t){!all(is.na(t))})]}
  return(as.data.frame(out))
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

getTF.gene_symbol <- function(TF.motif, HOCOMOCO.metadata, H.version="H12"){
  if (H.version=="H12"){gene_symbol <- filter(HOCOMOCO.metadata, name %in% TF.motif) %>% pull(masterlist_info.species.MOUSE.gene_symbol)} else {gene_symbol <- filter(HOCOMOCO.metadata, Model %in% TF.motif) %>% pull(`Transcription factor`)}
  return(gene_symbol)
}


plotTFfootprint.heatmap <- function(footprint.matrix, filter.unbound=FALSE, with.expression=FALSE, filter.no.expressed.TF=FALSE, Seurat.dataset=NULL, HOCOMOCO.metadata=NULL, H.version="H12"){
  require(ComplexHeatmap)
  require(circlize)
  require(patchwork)
  
  conditions <- str_remove(string = rownames(footprint.matrix), pattern = "_bound") %>% str_remove(pattern = "_score")
  rownames(footprint.matrix) <- conditions
  if (H.version=="H12"){colnames(footprint.matrix) <- str_remove(string=colnames(footprint.matrix), pattern = "_.*")} else {colnames(footprint.matrix) <- str_remove(string = colnames(footprint.matrix), pattern = "^.*\\.[:alpha:]+_")}
  if (filter.unbound){keep.i <- which(colSums(footprint.matrix)>0); footprint.matrix<-footprint.matrix[,keep.i]}

  if (with.expression & !is.null(Seurat.dataset) & !is.null(HOCOMOCO.metadata)){
    TF.motif.names <- colnames(footprint.matrix)
    TF.gene_symbols <- getTF.gene_symbol(TF.motif = TF.motif.names, HOCOMOCO.metadata = HOCOMOCO.metadata, H.version = H.version)
    TF.ensg <- convert_feature_identity(Seurat.dataset, assay = "RNA", features = TF.gene_symbols, feature.format = "symbol")
 
    DefaultAssay(Seurat.dataset) <- "RNA"
    TF.expression.data.avg <- AverageExpression(Seurat.dataset, assays = "RNA", features = na.omit(TF.ensg), group.by = "rv2.lineage")[[1]][,conditions]
    rownames(TF.expression.data.avg) <- convert_feature_identity(Seurat.dataset, assay = "RNA", features = rownames(TF.expression.data.avg), feature.format = "ens")
    
    TF.i <- match(rownames(TF.expression.data.avg), TF.gene_symbols)
    
    footprint.matrix <- footprint.matrix[,TF.i]
    
    # Filter non-expressed TFs
    
    if (filter.no.expressed.TF){
      TF.i.2 <- rowMeans(TF.expression.data.avg)>quantile(TF.expression.data.avg, 0.5)
      footprint.matrix <- footprint.matrix[,TF.i.2]
      TF.expression.data.avg<-TF.expression.data.avg[TF.i.2,]
    }
    TF.expression.data.avg<-t(scale(t(TF.expression.data.avg), center = FALSE))
    color.func.exp <- colorRamp2(c(min(TF.expression.data.avg), mean(apply(TF.expression.data.avg,1,mean)), max(TF.expression.data.avg)), c("white", "steelblue", "red"))
    
    # Perform TF expression splitting
    # TODO: Add iterative test to find number of centers?
    #kmeans.centers <- matrix(nrow = 4, ncol = ncol(TF.expression.data.avg),byrow = TRUE,c(seq(max(TF.expression.data.avg),0,length.out=ncol(TF.expression.data.avg)),seq(0,max(TF.expression.data.avg),length.out=ncol(TF.expression.data.avg)),rep(0,ncol(TF.expression.data.avg)),rep(mean(TF.expression.data.avg),ncol(TF.expression.data.avg))))
    TF.expression.mat.split <- kmeans(TF.expression.data.avg,centers=5, iter.max = 100)
  }
  
  #foot.dist <- dist(t(footprint.matrix), method="euclidean")
  #foot.hclust <- hclust(foot.dist, method = "complete")
  
  if (!any(footprint.matrix>0 & footprint.matrix<1)){color.func.foot <- c("black","green")} else {color.func.foot <- colorRamp2(c(min(footprint.matrix), mean(apply(footprint.matrix,2,mean)), max(footprint.matrix)), c("white", "yellow", "red"))}
  p.foot <- Heatmap(t(footprint.matrix), cluster_columns = FALSE, cluster_rows = TRUE, row_names_gp = gpar(fontsize = 6), col = color.func.foot, split=TF.expression.mat.split$cluster)
  
  if (with.expression){
    p.exp <- Heatmap(TF.expression.data.avg, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize = 6), col = color.func.exp,split=TF.expression.mat.split$cluster)
    p.final <- p.foot + p.exp + plot_layout(ncol=2)
    return(p.final)
  } else {
    return(p.foot)
  }
    
  
}