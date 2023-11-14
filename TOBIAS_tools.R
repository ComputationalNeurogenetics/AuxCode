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
  if(any(length(group)>1 | length(TF)>1)){
    errorCondition("Only one group and TF at the time are supported")
  } else {
    require(valr)
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

get.bounds <- function(tobias.gr, conditions, TF.motif, gr.filter=NULL){
  require(Signac)
  tob.sub <- tobias.gr[[TF.motif]]
  if (!is.null(gr.filter)){
    if (!class(gr.filter)=="GRanges"){gr.filter <- StringToGRanges(gr.filter)}
    tob.sub <- tob.sub[tob.sub %over% gr.filter]
  }
  bound.info <- elementMetadata(tob.sub)[,paste(conditions,"_bound",sep="")]
  out <- cbind(coordinate=GRangesToString(tob.sub), bound.info)
  return(out)
}

formFootprintMatrix.overConditions <- function(tobias.gr, conditions, gr.filter){
  all.motifs <- get.TF.motifs(tobias.gr)
  
  sapply(all.motifs, function(m){})
}
