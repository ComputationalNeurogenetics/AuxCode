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

get.bounds <- function(tobias.gr, conditions, TF.motif){
  require(Signac)
  tob.sub <- tobias.gr[[TF.motif]]
  bound.info <- elementMetadata(tob.sub)[,paste(conditions,"_bound",sep="")]
  out <- cbind(tibble(coordinate=GRangesToString(tob.sub)), bound.info)
  return(out)
}


TOBIAS.list.to.tibble <- function(tobias.list){
  # Assuming identical conditions for each element
  all.tf.motifs <- names(tobias.list)
  metadata.colnames <- colnames(elementMetadata(tobias.list[[1]]))
  conditions <- str_extract(string=metadata.colnames, pattern=".*_bound") %>% str_remove(pattern = "_bound") %>% na.omit()
  
}