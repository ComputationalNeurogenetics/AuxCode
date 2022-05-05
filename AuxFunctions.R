# Some additional functions ----

scale.by.acc <- function(tf.idf.mat, acc.mat, dummy=FALSE){
  if (dummy){
    return(tf.idf.mat)
  }
  feat.ids <- str_extract(string = colnames(tf.idf.mat), pattern = "chr[[:digit:]]{1,2}-[[:digit:]]+-[[:digit:]]+")
  cluster.ids <- str_extract(string = rownames(tf.idf.mat), pattern = "[[:digit:]]{1,2}")

  clust.ids.i <- match(cluster.ids, as.numeric(acc.mat$ident))
  tmp.acc <- data.frame(rep(acc.mat[,unique(feat.ids)], times=as.numeric(table(feat.ids))))[clust.ids.i,]
  tf.idf.mat*as.matrix(tmp.acc)
}

tf.idf.matrix <-function(data.matrix,log.tf=FALSE){
  idf <- apply(data.matrix,2,function(x){log(length(x)/sum(x>0))})
  if (log.tf){
    tf.idf.mat <- t(apply(data.matrix,1,function(x){((log(1+x)))*idf}))
  } else {
    tf.idf.mat <- t(apply(data.matrix,1,function(x){((x/sum(x)))*idf}))
  }
  return(tf.idf.mat)
}

TOBIAS.heatmap.plotter <- function(s.data, TFBS.data, genes, conditions, TF.meta.data, range.width, TOBIAS.format="SNAKEMAKE", TF.meta.format="HOCOMOCO", TF.filt=NULL){
  # Wrapper function to plot multiple TOBIAS heatmaps over conditions and genes from one dataset
  cond <- conditions
  cond.number <- as.numeric(str_remove(string=cond, pattern="^.*\\."))
  features.gr <- StringToGRanges(rownames(s.data))
  # Looping over genes
  for (g in genes){
    print(paste("Processing gene ", g, " of genes ", paste(genes, collapse = ","), sep=""))
    gene.coords <- range(Annotation(s.data)[Annotation(s.data)$gene_name==g], ignore.strand=TRUE)
    range.start <- start(gene.coords)
    range.end <- end(gene.coords)
    range.width <- range.width
    chr <- as.character(seqnames(gene.coords)@values) %>% str_remove(pattern = "[:lower:]{3}")

    gene.region <- construct_range(chr,range.start,range.end,range.width)
    features.in.region <- GRangesToString(features.gr[features.gr %over% StringToGRanges(gene.region)])

    # Looping over conditions
    for (cond in conditions){
      print(paste("Processing condition ", cond, " of conditions ", paste(conditions, collapse=","), sep=""))
      cond.number <- as.numeric(str_remove(string=cond, pattern="^.*\\."))
      TF.motifs.per.feat <- TF.motifs.per.feature.snakemake(features=features.in.region, TFBS.data=TFBS.data, features.in.region=features.in.region, region=gene.region, min.footprint.score=NULL, condition=cond)

      TF.motifs.per.feat$acc <- as.matrix(FetchData(s.data, vars=features.in.region, cells = WhichCells(s.data, idents=cond.number)))

      s.data.subset <- subset(s.data, cells=WhichCells(s.data, idents=cond.number))
      TF.motifs.per.feat$expr.pres <- find.TF.expr(TF.motifs.per.feat, s.data=s.data.subset, TF.metadata=TF.meta.data, TF.meta.format = TF.meta.format, TOBIAS.format = TOBIAS.format)

      s.data.subset <- tryCatch({
        LinkPeaks(object = s.data.subset,
                  peak.assay="peaks",
                  expression.assay="RNA_name",
                  expression.slot = "data",
                  gene.coords = NULL,
                  distance = range.width,
                  min.distance = NULL,
                  min.cells = 10,
                  method = "pearson",
                  genes.use = g,
                  n_sample = 200,
                  pvalue_cutoff = 0.05,
                  score_cutoff = 0.05,
                  verbose = TRUE
        )},
          error = function(e){
            s.data.subset
        }
      )

      TF.plot <- TF.heatmap(TF.mat.1 = TF.motifs.per.feat, TF.filt = TF.filt, TF.families = NULL, cluster.names = c("Footprint score"),expr.cutoff=0.1,TF.exprs=TRUE, row.cluster = TRUE, links.data=list(Links(s.data.subset)))
      tryCatch({
        draw(TF.plot, column_title = paste("E14, region.cluster: ", cond," ", g," region features with TF binding events",sep=""),column_title_gp = gpar(fontsize = 24))
      }, error=function(e){
        print("Skipping drawing as no bindings detected")
      })
    }


    }
}

find.TF.expr <- function(TF.footprint.data, s.data, TF.metadata, TF.meta.format, TOBIAS.format="TOBIAS"){
  DefaultAssay(s.data) <- "RNA_name"
  gene.name.space <- rownames(s.data)
  TF.ids <- rownames(TF.footprint.data$per.feat.mat)
  # Iterate through all TF motif ids in TF footprint data
  TF.gene.exprs <- sapply(TF.ids, function(motif.id){

    if (TF.meta.format=="JASPAR"){
      # Look for official gene name for the motif from JASPAR db derived TF metadata
      gene.name.tmp <- str_to_title(str_to_lower(filter(TF.metadata, motif.ids==motif.id) %>% pull(gene.name)))
    } else if (TF.meta.format=="HOCOMOCO") {
      # Look for official gene name for the motif from Hocomoco metadata derived via Motif Db
      if (TOBIAS.format=="TOBIAS"){
        gene.name.tmp <- str_to_title(filter(TF.metadata, motif==str_remove(string = motif.id, pattern = "^_")) %>% pull(geneSymbol))
      } else if (TOBIAS.format=="SNAKEMAKE"){
        gene.name.tmp <- str_to_title(filter(TF.metadata, Model==str_remove(string = motif.id, pattern = "^.*\\.[:upper:]_")) %>% pull(geneSymbol))
      }
    }

    if (gene.name.tmp %in% gene.name.space) {
      # If gene name present in gene.name.space we return its mean expression in the cell group within s.data
      gene.expr <- colMeans(FetchData(s.data, vars=gene.name.tmp))
      return(gene.expr)
    } else if (grepl(x=gene.name.tmp, pattern=".*::.*")) {
      # If gene name contains :: then there are two TFs binding one motif together, in that case we calculate mean expression for both genes from data in s.data, assuming they exist in the name space
      gene.names.tmp <- unlist(str_split(string=gene.name.tmp, pattern="::"))
      gene.names.tmp <- gene.names.tmp[gene.names.tmp %in% gene.name.space]
      gene.expr <- colMeans(FetchData(s.data, vars=gene.names.tmp))

      if (length(gene.expr)>1) {
        # If there are mean expression values for more than one gene we check if either is zero, then both are set to zero, otherwise we take mean of those.
        if (any(gene.expr==0)) {
          return(0)
        } else {
          return(mean(gene.expr))
        }
      } else {
        # If only one TF name was present in the gene name space then just return its mean value
        return(gene.expr)
      }
    } else {
      # If TF does not appear in gene name space we will return NA as its mean expression value
      return(NA)
    }
  })
  names(TF.gene.exprs) <- TF.ids
  return(TF.gene.exprs)
}

TF.heatmap <- function(TF.mat.1=NULL, TF.families=NULL, TF.filt=NULL, cluster.names=NA, links.data=NULL, TF.exprs=FALSE, expr.cutoff=NULL, row.cluster=FALSE, clustering_distance_rows="euclidean",clustering_method_rows = "complete"){
    TF.used <- rownames(TF.mat.1$per.feat.mat)
    if (!is.null(TF.filt)){
      TF.used <- TF.used[TF.used %in% paste(TF.filt,TF.filt,sep="_")]
    }
    TF.mat.1.expr <- TF.mat.1$expr.pres[TF.used]

    if (!is.null(expr.cutoff)){
      TF.expr.filt.l <- (TF.mat.1.expr > expr.cutoff) & !is.na(TF.mat.1.expr)
      TF.expr.filt.names <- names(TF.expr.filt.l)[TF.expr.filt.l==TRUE]
      TF.mat.1.expr <- TF.mat.1.expr[TF.expr.filt.l]
      TF.used <- TF.used[TF.used %in% names(TF.mat.1.expr)]
    }

    TF.mat.to.plot <- TF.mat.1$per.feat.mat[TF.used,]

    if (!max(TF.mat.to.plot)==0){

      col_fun = colorRamp2(c(0, max(TF.mat.to.plot)), c("white", "darkgreen"))

      if (!is.null(TF.families) & row.cluster==FALSE){
        row.split <- TF.families[rownames(TF.mat.to.plot)]
      } else {
        row.split <- NULL
      }

      # Format TF expression data into RowAnnotation if TF.exprs is TRUE
      if (TF.exprs){
        row_ha <- rowAnnotation(expr = anno_barplot(TF.mat.1$expr.pres[TF.used]))
        row_ha <- re_size(row_ha,width=unit(1,"inch"))
      } else {
        row_ha <- NULL
      }

      # Format Links data to col_ha if present and overlap gene region
      if (!is.null(links.data) & length(links.data[[1]])>0){
        overlapping.links <- any(StringToGRanges(links.data[[1]]$peak) %over% StringToGRanges(colnames(TF.mat.to.plot))==TRUE)
        if (overlapping.links){
          scores <- rep(0, ncol(TF.mat.1$acc))
          names(scores) <- colnames(TF.mat.1$acc)

          scores.tmp <- links.data[[1]]$score
          names(scores.tmp) <- links.data[[1]]$peak

          scores.tmp <- scores.tmp[names(scores.tmp) %in% names(scores)]

          scores[names(scores.tmp)] <- scores.tmp
          col_ha <- columnAnnotation(acc=anno_boxplot(TF.mat.1$acc, height = unit(4, "cm")), links=anno_barplot(scores, height = unit(4, "cm")))
        }
      } else {
        col_ha <- columnAnnotation(acc=anno_boxplot(TF.mat.1$acc, height = unit(4, "cm")))
      }

      TF.1.plot <- Heatmap(TF.mat.to.plot, cluster_rows = row.cluster, cluster_columns = FALSE, show_row_dend = TRUE, row_names_gp = gpar(fontsize = 6), col=col_fun, row_split=row.split, border = TRUE, row_title_rot = 0, row_gap = unit(2, "mm"), column_names_side = "top", heatmap_legend_param=list(title=cluster.names[1]), bottom_annotation = col_ha,  right_annotation = row_ha, clustering_distance_rows=clustering_distance_rows, clustering_method_rows=clustering_method_rows)
      return(TF.1.plot)
    } else {
      print("No bindings detected")
    }

}


TF.heatmap.diff <- function(TF.mat.1=NULL, TF.mat.2=NULL, TF.families=NULL, cluster.names=NA, links.data=NULL, TF.exprs=FALSE, expr.cutoff=NULL){
  if (all(!is.null(c(TF.mat.1,TF.mat.2)))){
    # Draw differential plot
    TF.used.i <- find.combined.non.empty.i(TF.mat.1$per.feat.mat, TF.mat.2$per.feat.mat)
    TF.used <- names(TF.used.i[TF.used.i==TRUE])

    TF.mat.1.expr <- TF.mat.1$expr.pres[TF.used]
    TF.mat.2.expr <- TF.mat.2$expr.pres[TF.used]

    if (!is.null(expr.cutoff)){
      TF.expr.filt.l <- (TF.mat.1.expr > expr.cutoff | TF.mat.2.expr > expr.cutoff) & (!is.na(TF.mat.1.expr) & !is.na(TF.mat.2.expr))
      TF.expr.filt.names <- names(TF.expr.filt.l)[TF.expr.filt.l==TRUE]
      TF.mat.1.expr <- TF.mat.1.expr[TF.expr.filt.l]
      TF.mat.2.expr <- TF.mat.2.expr[TF.expr.filt.l]
      TF.used <- TF.used[TF.used %in% names(TF.mat.1.expr)]
    }

    max.per.feat.score <- max(TF.mat.1$per.feat.mat, TF.mat.2$per.feat.mat)
    if (max.per.feat.score==0){max.per.feat.score<-1}
    col_fun = colorRamp2(c(0, max.per.feat.score), c("white", "darkgreen"))

    # Plot 1
    TF.mat.to.plot <- TF.mat.1$per.feat.mat[TF.used,]

    if (!is.null(TF.families)){
      row.split <- TF.families[rownames(TF.mat.to.plot)]
    } else {
      row.split <- NULL
    }
    # Format Links data to col_ha if present
    if (!is.null(links.data) & length(links.data[[1]])>0){
      overlapping.links <- any(StringToGRanges(links.data[[1]]$peak) %over% StringToGRanges(colnames(TF.mat.to.plot))==TRUE)
      if (overlapping.links){
        scores <- rep(0, ncol(TF.mat.1$acc))
        names(scores) <- colnames(TF.mat.1$acc)

        scores.tmp <- links.data[[1]]$score
        names(scores.tmp) <- links.data[[1]]$peak

        scores.tmp <- scores.tmp[names(scores.tmp) %in% names(scores)]

        scores[names(scores.tmp)] <- scores.tmp
        # TODO: This seems to fail if all scores are 0, needs to be handled properly
        col_ha <- columnAnnotation(acc=anno_boxplot(TF.mat.1$acc, height = unit(4, "cm")), links=anno_barplot(scores, height = unit(4, "cm")))
        }
      } else {
        col_ha <- columnAnnotation(acc=anno_boxplot(TF.mat.1$acc, height = unit(4, "cm")))
    }

    # Format TF expression data into RowAnnotation if TF.exprs is TRUE
    if (TF.exprs){
      #browser()
      row_ha <- rowAnnotation(expr = anno_barplot(TF.mat.1.expr))
    } else
      {
      row_ha <- NULL
    }

    TF.1.plot <- Heatmap(TF.mat.to.plot, cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col=col_fun, row_split=row.split, border = TRUE, row_title_rot = 0, row_gap = unit(2, "mm"), column_names_side = "top", column_title = cluster.names[1], heatmap_legend_param=list(title=cluster.names[1]), bottom_annotation = col_ha, right_annotation = row_ha)

    # Plot 2
    TF.mat.to.plot <- TF.mat.2$per.feat.mat[TF.used,]

    if (!is.null(TF.families)){
      row.split <- TF.families[rownames(TF.mat.to.plot)]
    } else {
      row.split <- NULL
    }
    # Format Links data to col_ha if present and overlap gene region
    if (!is.null(links.data) & length(links.data[[2]])>0){
      overlapping.links <- any(StringToGRanges(links.data[[2]]$peak) %over% StringToGRanges(colnames(TF.mat.to.plot))==TRUE)
      if (overlapping.links){
        scores <- rep(0, ncol(TF.mat.2$acc))
        names(scores) <- colnames(TF.mat.2$acc)

        scores.tmp <- links.data[[2]]$score
        names(scores.tmp) <- links.data[[2]]$peak

        scores.tmp <- scores.tmp[names(scores.tmp) %in% names(scores)]

        scores[names(scores.tmp)] <- scores.tmp
        col_ha <- columnAnnotation(acc=anno_boxplot(TF.mat.2$acc, height = unit(4, "cm")), links=anno_barplot(scores, height = unit(4, "cm")))
    }
    } else {
        col_ha <- columnAnnotation(acc=anno_boxplot(TF.mat.2$acc, height = unit(4, "cm")))
    }


    # Format TF expression data into RowAnnotation if TF.exprs is TRUE
    if (TF.exprs){
      row_ha <- rowAnnotation(expr = anno_barplot(TF.mat.2.expr))
    } else
    {
      row_ha <- NULL
    }

    TF.2.plot <- Heatmap(TF.mat.to.plot, cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col=col_fun, row_split=row.split, border = TRUE, row_title_rot = 0, row_gap = unit(2, "mm"), column_names_side = "top", column_title = cluster.names[2], heatmap_legend_param=list(title=cluster.names[2]), bottom_annotation = col_ha, right_annotation = row_ha)

    # Calculate differential
    TF.diff.mat <- TF.diff(TF.mat.1$per.feat.mat[TF.used,],TF.mat.2$per.feat.mat[TF.used,])

    TF.mat.to.plot <- TF.diff.mat
    col_fun = colorRamp2(c(min(TF.mat.to.plot), 0, max(TF.mat.to.plot)), c("blue", "white", "red"))

    if (!is.null(TF.families)){
      row.split <- TF.families[rownames(TF.mat.to.plot)]
    } else {
      row.split <- NULL
    }

    mean.acc.diff <- colMeans(TF.mat.1$acc)- colMeans(TF.mat.2$acc)
    col_ha <- columnAnnotation(acc.diff = anno_barplot(mean.acc.diff, height = unit(4, "cm"), gp = gpar(fill = ifelse(mean.acc.diff>0, "red", "blue"))))

    # Format TF expression data into RowAnnotation if TF.exprs is TRUE
    if (TF.exprs){
      #browser()
      diff.exp <- TF.mat.1$expr.pres[TF.used] - TF.mat.2$expr.pres[TF.used]
      row_ha <- rowAnnotation(expr = anno_barplot(diff.exp))
    } else
    {
      row_ha <- NULL
    }

    TF.3.plot <- Heatmap(TF.mat.to.plot, cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col=col_fun, row_split=row.split, border = TRUE, row_title_rot = 0, row_gap = unit(2, "mm"), column_names_side = "top",  column_title = "Difference", heatmap_legend_param=list(title="Diff (1st-2nd)"), bottom_annotation = col_ha, right_annotation = row_ha)
    TF.plot.combined <- TF.1.plot + TF.2.plot + TF.3.plot
    return(TF.plot.combined)
  } else {
    # Error
  }
}

construct_range <- function(chr,gene.start,gene.end, width){
  return(paste("chr",chr,"-",gene.start-width,"-",gene.end+width,sep=""))
}

find_TFBS_range <- function(tobias_set, region, filter.bound=F, return.empty=F){
  region.gr <- StringToGRanges(region)
  hits <- lapply(tobias_set,function(ts){

    # Look for tobias sets overlapping region of interest only if there is common sequence level, otherwise return empty granges set
    if (all(seqlevels(ts) %in% seqlevels(region.gr)==FALSE)){
      ts <- ts[FALSE,]
    } else {
    ts <- ts[ts %over% region.gr]
    }

    if (filter.bound & length(ts) > 0){
      tfbs.metadata.tmp <- colnames(ts@elementMetadata)
      bound.cols <- grep(tfbs.metadata.tmp, pattern = ".*bound")
      bound.i <- apply(ts@elementMetadata[,bound.cols],1,function(b){any(b==1)})
      return(ts[bound.i])
    } else {
      return(ts)
    }
  })

  ts.hits.count <- sapply(hits,length)
  if (return.empty){
    ts.hits <- hits
  } else {
    ts.hits <- hits[ts.hits.count>0]
  }
  return(ts.hits)
}

Find.TF.families <- function(TOBIAS_res){
  require(JASPAR2020)
  require(TFBSTools)
  TF.family <- unlist(sapply(names(TOBIAS_res), function(tf){
    opts <- list()
    opts[["ID"]] <- str_replace(string = tf, pattern = ".*_",replacement = "")
    TF.tmp.data <- getMatrixSet(JASPAR2020, opts)
    return(TF.tmp.data[[1]]@tags$family)
  }, simplify = TRUE))
  return(TF.family)
}

TF.diff <- function(TF.mat.1, TF.mat.2){
  TF.mat.diff <- TF.mat.1 - TF.mat.2
  return(TF.mat.diff)
}

find.combined.non.empty.i <- function(TF.matrix.1, TF.matrix.2){
  TF.1.i <- apply(TF.matrix.1,1,function(r){any(r>0)})
  TF.2.i <- apply(TF.matrix.2,1,function(r){any(r>0)})
  #browser()
  #joint.TF.motifs <- intersect(names(TF.1.i),names(TF.2.i))
  return(TF.1.i | TF.2.i)
}

TF.motifs.per.feature.snakemake <- function(features=NULL, TFBS.data, region, features.in.region=NULL, min.footprint.score=NULL, condition){
  # Define features dimension (cols)
  if (!is.null(features)){
    features.gr <- StringToGRanges(features)
    names(features.gr) <- rownames(features)
    features.in.region <- features.gr[features.gr %over% StringToGRanges(region)]
  } else {
    features.in.region <- StringToGRanges(features.in.region)
  }
  print(paste("Found ", length(features.in.region), " features in the region", sep=""))
  if (!class(TFBS.data)=="CompressedGRangesList"){
    TFBS.data <- GRangesList(TFBS.data)
  }

  TFBS.in.features <- lapply(TFBS.data, function(tfbs){

    # Subset granges based on bound==1 on given condition
    tfbs.colnames <- colnames(mcols(tfbs))
    tfbs.i <- which(tfbs.colnames==paste(condition,"_bound",sep=""))
    tfbs.filt <- tfbs[elementMetadata(tfbs)[,tfbs.i]==1,]
    tmp.hits <- findOverlaps(query = tfbs.filt, subject = features.in.region, minoverlap = 1)
    tmp.features <- GRangesToString(features.in.region[subjectHits(tmp.hits)])
    tfbs.hits <- tfbs[queryHits(tmp.hits)]
    if (length(tfbs.hits)>0){
      tfbs.hits$feature <- tmp.features
    }
    return(tfbs.hits)
  })

  TF.hit.count <- sapply(TFBS.in.features, length)
  TF.hits <- TFBS.in.features[TF.hit.count>0]

  print(paste("Found ", sum(TF.hit.count>0), " hits in the region", sep=""))

  # Create zero matrix
  TF.motif.matrix <- matrix(0, nrow = length(TFBS.data), ncol=length(features.in.region))
  rownames(TF.motif.matrix) <- names(TFBS.data)
  colnames(TF.motif.matrix) <- GRangesToString(features.in.region)

  # Loop over all TFBS binding events which overlapped features in the gene region
  lapply(names(TF.hits), function(tf){
    feature <- TF.hits[[tf]]$feature
    tfbs.colnames <- colnames(mcols(TF.hits[[tf]]))
    tfbs.i <- which(tfbs.colnames==paste(condition,"_score",sep=""))
    footprint.scores <- TF.hits[[tf]][,tfbs.i]

    avg.footprint.score.per.feat <- tapply(INDEX=feature, X=mcols(footprint.scores)[,1], FUN=mean)
    TF.motif.matrix[tf,names(avg.footprint.score.per.feat)] <<- avg.footprint.score.per.feat
  })
  return(list(per.feat.mat=TF.motif.matrix))
}


TF.motifs.per.feature <- function(features, TFBS.data, region, min.footprint.score=NULL){
  # Define features dimension (cols)
  features.gr <- StringToGRanges(features)
  names(features.gr) <- rownames(features)
  features.in.region <- features.gr[features.gr %over% StringToGRanges(region)]
  print(paste("Found ", length(features.in.region), " features in the region", sep=""))

  TFBS.gr.list <- GRangesList(TFBS.data)
  TFBS.in.features <- lapply(TFBS.gr.list, function(tfbs){
    tmp.hits <- findOverlaps(query = tfbs, subject = features.in.region, minoverlap = 1)
    tmp.features <- GRangesToString(features.in.region[subjectHits(tmp.hits)])
    tfbs.hits <- tfbs[queryHits(tmp.hits)]
    if (length(tfbs.hits)>0){
      tfbs.hits$feature <- tmp.features
    }
    return(tfbs.hits)
  })

  TF.hit.count <- sapply(TFBS.in.features, length)
  TF.hits <- TFBS.in.features[TF.hit.count>0]

  # TODO: Add print for found tfbs in features

  # Create zero matrix
  TF.motif.matrix <- matrix(0, nrow = length(TFBS.data), ncol=length(features.in.region))
  rownames(TF.motif.matrix) <- names(TFBS.data)
  colnames(TF.motif.matrix) <- GRangesToString(features.in.region)

  # Loop over all TFBS binding events which overlapped features in the gene region
  lapply(names(TF.hits), function(tf){
    feature <- TF.hits[[tf]]$feature
    avg.footprint.score.per.feat <- tapply(INDEX=feature, X=TF.hits[[tf]]$footprint_score, FUN=mean)
    TF.motif.matrix[tf,names(avg.footprint.score.per.feat)] <<- avg.footprint.score.per.feat
  })
  return(list(per.feat.mat=TF.motif.matrix))
}


get_BINDetect_snakemake_results <- function(res_path,parallel=F, mc.cores=NULL){

  #'@param res_path (str): Path to the folder where TOBIAS BINDetect results are stored
  #'
  #'@returns a named list (Large list) consisting of granges for each result sub folder in @param res_path.
  #'         The list can be conveniently accessed, for example, with out_list$gene_TFBSname
  #'
  #'@example get_BINDetect_snakemake_results("/path/to/TOBIAS_framework/outs/TFBS/")
  #'
  # 'Dependencies'
  library(GenomicRanges)
  library(magrittr)

  # Reject the .txt, .pdf, etc. files with regex.
  # Apparently all sub folders are of form gene_TFBSname.n where n \in {1,2,3}
  motif.res.folders <- list.files(res_path, pattern = "(.*\\.H[0-9]{2}MO\\.[A-Z]{1})|(\\.[0-9])")

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
    # TODO: This could be significantly faster if one first catenates all files and reads it as one pass into df and then GR
    GenomicRanges::makeGRangesFromDataFrame(overview.df, keep.extra.columns = TRUE,
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
      GenomicRanges::makeGRangesFromDataFrame(overview.df, keep.extra.columns = TRUE,
                                              seqnames.field = "TFBS_chr",
                                              start.field = "TFBS_start",
                                              end.field = "TFBS_end",
                                              strand.field = "TFBS_strand")
    }, mc.cores=mc.cores)
    names(out_list) <- motif.res.folders


  }
  return(out_list)
}

get_BINDetect_results <- function(res_path, col.names="Default") {
  #'
  #' Queries TOBIAS BINDetect result folder and selects result bed files.
  #'
  #' -------------------------Briefly in pseudocode------------------------------
  #'
  #' Access the BINDetect result folder
  #'
  #' Initialize a result list
  #'
  #' For each sub folder in the result folder do:
  #'     Open the sub folder (henceforth sub folder = gene_TFBSname)
  #'     Read the gene_TFBSname/beds/gene_TFBS_name_bound.bed
  #'     Convert the bed file into a granges
  #'     Append the granges to the result list
  #'
  #' Name the list after string vector (gene_TFBSname1, gene_TFBSname2,...)
  #'
  #' Return the result list
  #' ----------------------------------------------------------------------------
  #'
  #' Requires installation of packages:
  #'     GenomicRanges
  #'     Magrittr
  #'
  #'@param res_path (str): Path to the folder where TOBIAS BINDetect results are stored
  #'
  #'@returns a named list (Large list) consisting of granges for each result sub folder in @param res_path.
  #'         The list can be conveniently accessed, for example, with out_list$gene_TFBSname
  #'
  #'@example get_BINDetect_results("/path/to/TOBIAS_framework/outs/TOBIAS_BINDetect_output/")
  #'
  # 'Dependencies'
  library(GenomicRanges)
  library(magrittr)

  # Reject the .txt, .pdf, etc. files with regex.
  # Apparently all sub folders are of form gene_TFBSname.n where n \in {1,2,3}
  filenames <- list.files(res_path, pattern = "(.*\\.H[0-9]{2}MO\\.[A-Z]{1})|(\\.[0-9])")
  # The actual loop as described in pseudo
  out_list <- lapply(filenames, function(name) {
    # Access the sub folder's contents.
    # This should be of form res_path/gene_TFBSname.n/beds/
    bound.bed.path <- paste0(res_path, name) %>% paste0("/beds/")
    # Assuming that the files are in constant order where for index set i = {1,2,3} in the folder ./beds/:
    #   1.   gene_TFBSname.n_all.bed
    #   2.   gene_TFBSname.n_bound.bed
    #   3.   gene_TFBSname.n_unbound.bed
    bound.bed.file <- list.files(bound.bed.path)[2]

    # Merge the former two to get full path to the bed
    bound.bed.full.path <- paste0(bound.bed.path, bound.bed.file)

    # A little derail, but apparently the most simple way to name each column in the granges is
    # to convert the bed-file into a column-named data frame.
    # The Granges inherits the column names and thus is can be indexed by column names.
    bound.bed.df <- data.frame(read.table(bound.bed.full.path))
    # Df column names created after column names in respath/gene_TFBSname.n/gene_TFBSname.n_overview.txt
    if (col.names[1]=="Default"){
    col.names.in.input <- c("TFBS_chr",    "TFBS_start", "TFBS_end",
                            "TFBS_name","TFBS_score", "TFBS_strand",
                            "peak_chr", "peak_start",   "peak_end",
                            "peak_id", "peak_score", "peak_strand",
                            "gene_id", "footprint_score")
    colnames(bound.bed.df) <- col.names.in.input
    } else {
      col.names.in.input <- col.names
      colnames(bound.bed.df) <- col.names.in.input
    }
    # Conversion into a Granges object. The keep.extra.columns argument stores
    # other columns into the metadata slot, by name.
    bound.bed.granges <- makeGRangesFromDataFrame(bound.bed.df,
                                                  keep.extra.columns = TRUE,
                                                  seqnames.field = "TFBS_chr",
                                                  start.field = "TFBS_start",
                                                  end.field = "TFBS_end",
                                                  strand.field = "TFBS_strand")

    # Return the list of Granges objects
    return(bound.bed.granges)
  })
  # Name each Granges with corresponding gene_TFBSname.n
  names(out_list) <- list.files(res_path, pattern = "(.*\\.H[0-9]{2}MO\\.[A-Z]{1})|(\\.[0-9])")

  # Return the Granges list
  return(out_list)

}

get_TFBS_overview_results <- function(res_path) {
  #'
  #' Queries TOBIAS BINDetect result folder and finds TFBS overview of each TF.
  #'
  #' -------------------------Briefly in pseudocode------------------------------
  #'
  #' Access the BINDetect result folder
  #'
  #' Initialize a result list
  #'
  #' For each sub folder in the result folder do:
  #'     Open the sub folder (henceforth sub folder = gene_TFBSname)
  #'     Read the gene_TFBSname/gene_TFBS_name_overview.txt
  #'     Convert the txt file into a granges
  #'     Append the granges to the result list
  #'     Convert list to Grangeslist
  #'
  #' Name the list after string vector (gene_TFBSname1, gene_TFBSname2,...)
  #'
  #' Return the result list
  #' ----------------------------------------------------------------------------
  #'
  #' Requires installation of packages:
  #'     GenomicRanges
  #'     Magrittr
  #'
  #'@param res_path (str): Path to the folder where TOBIAS BINDetect results are stored
  #'
  #'@returns a named list (Large list) consisting of granges for each result sub folder in @param res_path.
  #'         The list can be conveniently accessed, for example, with out_list$gene_TFBSname
  #'
  #'@example get_BINDetect_results("/path/to/TOBIAS_framework/outs/TOBIAS_BINDetect_output/")
  #'
  # 'Dependencies'
  library(GenomicRanges)
  library(magrittr)

  # Reject the .txt, .pdf, etc. files with regex.
  # Apparently all sub folders are of form gene_TFBSname.n where n \in {1,2,3}
  filenames <- list.files(res_path, pattern = "\\.[0-9]")
  # The actual loop as described in pseudo
  out_list <- lapply(filenames, function(name) {
    # Access the sub folder's contents.
    # This should be of form res_path/gene_TFBSname.n/beds/
    tfbs.overview.path <- paste0(res_path, name,"/") %>% paste0(name,"_overview.txt",sep="")

    # A little derail, but apparently the most simple way to name each column in the granges is
    # to convert the bed-file into a column-named data frame.
    # The Granges inherits the column names and thus is can be indexed by column names.
    tfbs.overview.df <- data.frame(read_tsv(tfbs.overview.path, col_names = TRUE))

    # Conversion into a Granges object. The keep.extra.columns argument stores
    # other columns into the metadata slot, by name.
    tfbs.overview.granges <- makeGRangesFromDataFrame(tfbs.overview.df,
                                                  keep.extra.columns = TRUE,
                                                  seqnames.field = "TFBS_chr",
                                                  start.field = "TFBS_start",
                                                  end.field = "TFBS_end",
                                                  strand.field = "TFBS_strand")

    # Return the list of Granges objects
    return(tfbs.overview.granges)
  })
  # Name each Granges with corresponding gene_TFBSname.n
  names(out_list) <- list.files(res_path, pattern = "\\.[0-9]")
  # Return the Granges list
  return(out_list)

}


RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  obj[['RNA_name']] <- obj[['RNA']]
  RNA <- obj@assays$RNA_name

  #tmp.conv <- tibble(id=RNA@counts@Dimnames[[1]], symbol=newnames)

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    #if (length(RNA@scale.data)) dimnames(RNA@scale.data)[[1]]    <- tmp.conv$symbol[match(dimnames(RNA@scale.data)[[1]],tmp.conv$id)]
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA_name <- RNA
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)

#
# LinkPeaksId <- function (object, peak.assay, expression.assay, expression.slot = "data",
#                          gene.coords = NULL, distance = 5e+05, min.distance = NULL,
#                          min.cells = 10, method = "pearson", genes.use = NULL, n_sample = 200,
#                          pvalue_cutoff = 0.05, score_cutoff = 0.05, verbose = TRUE)
# {
#   if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
#     stop("The requested assay is not a ChromatinAssay")
#   }
#   if (!is.null(x = min.distance)) {
#     if (!is.numeric(x = min.distance)) {
#       stop("min.distance should be a numeric value")
#     }
#     if (min.distance < 0) {
#       warning("Requested a negative min.distance value, setting min.distance to zero")
#       min.distance <- NULL
#     }
#     else if (min.distance == 0) {
#       min.distance <- NULL
#     }
#   }
#   if (is.null(x = gene.coords)) {
#     gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = object[[peak.assay]]))
#   }
#   meta.features <- GetAssayData(object = object, assay = peak.assay,
#                                 slot = "meta.features")
#   features.match <- c("GC.percent", "count")
#   if (!("GC.percent" %in% colnames(x = meta.features))) {
#     stop("GC content per peak has not been computed.\n",
#          "Run RegionsStats before calling this function.")
#   }
#   peak.data <- GetAssayData(object = object, assay = peak.assay,
#                             slot = "counts")
#   if (!("count" %in% colnames(x = meta.features))) {
#     hvf.info <- FindTopFeatures(object = peak.data)
#     hvf.info <- hvf.info[rownames(x = meta.features), ]
#     meta.features <- cbind(meta.features, hvf.info)
#   }
#   expression.data <- GetAssayData(object = object, assay = expression.assay,
#                                   slot = expression.slot)
#   peakcounts <- meta.features[rownames(x = peak.data), "count"]
#   genecounts <- rowSums(x = expression.data > 0)
#   peaks.keep <- peakcounts > min.cells
#   genes.keep <- genecounts > min.cells
#   peak.data <- peak.data[peaks.keep, ]
#   if (is.null(x = genes.use)) {
#     expression.data <- expression.data[genes.keep, ]
#   }
#   else {
#     genes.keep <- intersect(x = names(x = genes.keep[genes.keep]),
#                             y = genes.use)
#     expression.data <- expression.data[genes.keep, , drop = FALSE]
#   }
#   if (verbose) {
#     message("Testing ", nrow(x = expression.data), " genes and ",
#             sum(peaks.keep), " peaks")
#   }
#   genes <- rownames(x = expression.data)
#   gene.coords.use <- gene.coords[gene.coords$gene_id %in%
#                                    genes, ]
#   peaks <- granges(x = object[[peak.assay]])
#   peaks <- peaks[peaks.keep]
#   peak_distance_matrix <- DistanceToTSS(peaks = peaks, genes = gene.coords.use,
#                                         distance = distance)
#   if (!is.null(x = min.distance)) {
#     peak_distance_matrix_min <- DistanceToTSS(peaks = peaks,
#                                               genes = gene.coords.use, distance = min.distance)
#     peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
#   }
#
#   colnames(peak_distance_matrix) <- genes.use
#   if (sum(peak_distance_matrix) == 0) {
#     stop("No peaks fall within distance threshold\n", "Have you set the proper genome and seqlevelsStyle for ",
#          peak.assay, " assay?")
#   }
#   genes.use <- colnames(x = peak_distance_matrix)
#   all.peaks <- rownames(x = peak.data)
#   peak.data <- t(x = peak.data)
#   coef.vec <- c()
#   gene.vec <- c()
#   zscore.vec <- c()
#   if (nbrOfWorkers() > 1) {
#     mylapply <- future_lapply
#   }
#   else {
#     mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
#   }
#   res <- mylapply(X = seq_along(along.with = genes.use), FUN = function(i) {
#     peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
#     gene.expression <- t(x = expression.data[genes.use[[i]],
#                                              , drop = FALSE])
#     gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
#     if (sum(peak.use) < 2) {
#       return(list(gene = NULL, coef = NULL, zscore = NULL))
#     }
#     else {
#       peak.access <- peak.data[, peak.use, drop = FALSE]
#       coef.result <- corSparse(X = peak.access, Y = gene.expression)
#       rownames(x = coef.result) <- colnames(x = peak.access)
#       coef.result <- coef.result[abs(x = coef.result) >
#                                    score_cutoff, , drop = FALSE]
#       if (nrow(x = coef.result) == 0) {
#         return(list(gene = NULL, coef = NULL, zscore = NULL))
#       }
#       else {
#         peaks.test <- rownames(x = coef.result)
#         trans.peaks <- all.peaks[!grepl(pattern = paste0("^",
#                                                          gene.chrom), x = all.peaks)]
#         meta.use <- meta.features[trans.peaks, ]
#         pk.use <- meta.features[peaks.test, ]
#         bg.peaks <- lapply(X = seq_len(length.out = nrow(x = pk.use)),
#                            FUN = function(x) {
#                              MatchRegionStats(meta.feature = meta.use,
#                                               query.feature = pk.use[x, , drop = FALSE],
#                                               features.match = c("GC.percent", "count",
#                                                                  "sequence.length"), n = n_sample, verbose = FALSE)
#                            })
#         bg.access <- peak.data[, unlist(x = bg.peaks),
#                                drop = FALSE]
#         bg.coef <- corSparse(X = bg.access, Y = gene.expression)
#         rownames(bg.coef) <- colnames(bg.access)
#         zscores <- vector(mode = "numeric", length = length(x = peaks.test))
#         for (j in seq_along(along.with = peaks.test)) {
#           coef.use <- bg.coef[(((j - 1) * n_sample) +
#                                  1):(j * n_sample), ]
#           z <- (coef.result[j] - mean(x = coef.use))/sd(x = coef.use)
#           zscores[[j]] <- z
#         }
#         names(x = coef.result) <- peaks.test
#         names(x = zscores) <- peaks.test
#         zscore.vec <- c(zscore.vec, zscores)
#         gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
#         coef.vec <- c(coef.vec, coef.result)
#       }
#       gc(verbose = FALSE)
#       pval.vec <- pnorm(q = -abs(x = zscore.vec))
#       links.keep <- pval.vec < pvalue_cutoff
#       if (sum(x = links.keep) == 0) {
#         return(list(gene = NULL, coef = NULL, zscore = NULL))
#       }
#       else {
#         gene.vec <- gene.vec[links.keep]
#         coef.vec <- coef.vec[links.keep]
#         zscore.vec <- zscore.vec[links.keep]
#         return(list(gene = gene.vec, coef = coef.vec,
#                     zscore = zscore.vec))
#       }
#     }
#   })
#   browser()
#   gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`,
#                                               1))
#   coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`,
#                                               2))
#   zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`,
#                                                 3))
#   if (length(x = coef.vec) == 0) {
#     if (verbose) {
#       message("No significant links found")
#     }
#     return(object)
#   }
#   peak.key <- seq_along(along.with = unique(x = names(x = coef.vec)))
#   names(x = peak.key) <- unique(x = names(x = coef.vec))
#   coef.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = coef.vec)],
#                               x = coef.vec, dims = c(length(x = genes.use), max(peak.key)))
#   rownames(x = coef.matrix) <- genes.use
#   colnames(x = coef.matrix) <- names(x = peak.key)
#   links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
#   z.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = zscore.vec)],
#                            x = zscore.vec, dims = c(length(x = genes.use), max(peak.key)))
#   rownames(x = z.matrix) <- genes.use
#   colnames(x = z.matrix) <- names(x = peak.key)
#   z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
#   links$zscore <- z.lnk$score
#   links$pvalue <- pnorm(q = -abs(x = links$zscore))
#   links <- links[links$pvalue < pvalue_cutoff]
#   Links(object = object[[peak.assay]]) <- links
#   return(object)
# }
#
# environment(LinkPeaksId) <- asNamespace("Signac")

convert_feature_identity <- function(object, assay, features, feature.format = "symbol") {

  #'
  #' Converts ENS ID -> gene symbol and vice versa
  #' Returns a vector of length(features) of either matches or NAs, in corresponding indices.
  #'
  #' Assumes libraries dplyr and seuratObject.
  #' Moreover, requires seuratObject[["assay"]] to contain df/tbl of ENS-symbol correspondences.
  #'
  # Protective tests
  if (!any(feature.format %in% c("ens", "symbol"))) {
    stop("Feature format should be a sting: either 'symbol' or 'ens'")
  }
  if (length(features) == 0) {
    stop("No features found in argument 'features'")
  }
  if (feature.format == "ens" & !all(grepl("*ENS", features))) {
    stop("Found non-ENS ID for argument feature format 'ens'")
  }
  if (feature.format == "symbol" & any(grepl("*ENS", features))) {
    stop("Found ENS ID for argument feature format 'symbol'")
  }
  # Diverging execution: case if provided features are ENSEBL IDs => conversion to symbols
  if (feature.format == "ens") {

    object.features <- object[[assay]][[]] %>%
      rownames_to_column(var = "gene_id") %>%
      as_tibble() %>%
      dplyr::select("gene_id", "feature_symbol") %>%
      dplyr::filter(gene_id %in% features)
    match.index <- match(features, object.features$gene_id, nomatch = NA)
    v.out <- sapply(match.index, function (i) { ifelse(is.na(i), NA, object.features$feature_symbol[i])})

    sprintf("Instance: Found matching for %d features out of total %d provided features", sum(!is.na(v.out)), length(features)) %>%
      print()

    return (v.out)
  }

  # Case: otherwise provided symbols => conversion to ENS IDs
  object.features <- object[[assay]][[]] %>%
    rownames_to_column(var = "gene_id") %>%
    as_tibble() %>%
    dplyr::select("gene_id", "feature_symbol") %>%
    dplyr::filter(feature_symbol %in% features)

  match.index <- match(features, object.features$feature_symbol, nomatch = NA)
  v.out <- sapply(match.index, function (i) { ifelse(is.na(i), NA, object.features$gene_id[i])})

  sprintf("Instance: Found matching for %d features out of total %d provided features", sum(!is.na(v.out)), length(features)) %>%
    print()

  return (v.out)
}

fit.GAM <- function(rna,p.time){
  # Fit GAM for each gene using pseudotime as independent variable.
  t <- p.time
  gam.pval <- apply(rna, 1, function(z){
    d <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    summary(tmp)[4][[1]][1,5]
  })
}

list2df_tibble <- function(x) {
  tmp <- purrr::map(x, tibble::as_tibble)
  dplyr::bind_rows(tmp, .id = "name")
}

# Function to transform separate txt files after MACS2 into R object to be passed then into Seurat/Signac. Note this is extremely
# IO heacy operation, run with fast storage
peaks_to_matrix <- function(sample.name = sample.name, path = path){
# Script to generate counts DF for from separate txt files for Seurat based pipeline downstream
# Recommended to be run somewhere with very fast filesystem

library(plyr)
library(data.table)

# Sample name
sample.name <- sample.name

#' Read new peak count data in
path <- path
files <- list.files(path,pattern = "\\.txt$")
length(files)

#' Assuming tab separated values with a header
datalist <- llply(files, function(x)fread(paste0(path,x))$V4, .progress='text')
datafr <- do.call("cbind", datalist)

# Separating cell barcodes from bam filenames to be used as observation (column) names. Not perhaps the most elegant solution but it works.
barcodes <- substr(gsub(pattern='^TAG_CB_', replacement='', x=sapply(strsplit(files,'\\.'),'[', 1)),start=0,stop=16)
colnames(datafr) <- barcodes

saveRDS(datafr, paste(sample.name,"_CountsDF.Rds",sep=""))
}

# This is from Staija lab vignette https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html
# topTFs <- function(celltype, padj.cutoff = 1e-2) {
#   ctmarkers_rna <- dplyr::filter(
#     markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>%
#     arrange(-RNA.auc)
#   ctmarkers_motif <- dplyr::filter(
#     markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>%
#     arrange(-motif.auc)
#   top_tfs <- inner_join(
#     x = ctmarkers_rna[, c(2, 11, 6, 7)],
#     y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
#   )
#   top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
#   top_tfs <- arrange(top_tfs, -avg_auc)
#   return(top_tfs)
# }


create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))))
}


smooth.TFz.sd.rand <- function(TFz, cell.order, f=2/3, rand.n=10){
  rand.sd<-sapply(1:rand.n, function(r){
    sd(lowess(y = TFz, x = sample(x = cell.order, size = length(cell.order)), f=f)$y)
  })
  return(rand.sd)
}



Find_TF_Marker_associations <- function(motifs, markers, db.name, adj.p.value = 0.05){
  # Checking requirements ----
  require(enrichR)
  require(tidyverse)

  # Calculate enrichments of a specific TF target genes among the markers ----
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE else websiteLive <- TRUE
  if (websiteLive) enriched <- enrichr(markers, db.name)
  enriched.tb <- as_tibble(enriched)
  enriched.tb.filt <- filter(enriched.tb[[1]],Adjusted.P.value <= adj.p.value)

  # Extract TF names from motifs ----
  TF.names.from.motifs <- gsub(pattern = "-MA.*", x = motifs,replacement = "")

  # Extract TF names from enrichments ----
  TF.names.from.enrichments <- gsub(pattern = "\\s+.*", x = pull(enriched.tb.filt,Term), replacement = "")

  # Intersect names ----
  intersected.TF.names <- intersect(TF.names.from.enrichments,TF.names.from.motifs)

  # Pull out genes present both in markers and enriched among the downstream target genes of each intersected TF ----
  program.genes <- lapply(intersected.TF.names, function(term){
    # Pulling out genes of each term from the intersect
    all.terms <- pull(enriched.tb.filt,Term)
    term.i <- grep(pattern = term, x = all.terms)
    genes.downstream <- strsplit(x = pull(enriched.tb.filt[term.i,], Genes),split=";")[[1]]
    # Intersecting genes downstream with markers
    genes.for.program <- intersect(genes.downstream,toupper(markers))
    # Add TF itself to the program
    str_to_title(c(genes.for.program, term))
  })
  # Format output ----
  names(program.genes) <- intersected.TF.names
  return(list(programs=program.genes, db.name=db.name, adj.p.value.threshold = adj.p.value))
}

Fetch_Enrichments <- function(genes,db.name, adj.p.value = 0.05){
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE else websiteLive <- TRUE
  if (websiteLive) {
    enriched <- enrichr(genes, db.name)
    enriched <- lapply(enriched, function(e){filter(as_tibble(e),Adjusted.P.value <= adj.p.value)})
  }
  enriched
}


search_replace <- function(v,replacement){
  corrected <- lapply(v, function(x){
    if (x %in% names(replacement)){
      replacement[[x]]
    } else {
      x
    }
  })
  return(unlist(corrected, use.names=FALSE))
}

motif.to.geneName <- function(motif, dash.type="-"){
  gene.part <- strsplit(motif,split=paste(dash.type,"MA",sep=""))[[1]][1]
  if (length(grep(pattern="::", gene.part)) == 0){
    # Strip (ver something) away
    gene.part <- gsub(pattern = "\\(.+\\)", x = gene.part, replacement = "")
    gene.part <- gsub(pattern = "var.\\d", x = gene.part, replacement = "")
    return(gene.part)
  } else {
    genes.part <- strsplit(gene.part,split = "::")[[1]]
     # Strip (ver something) away
    genes.part <- gsub(pattern = "\\(.+\\)", x = genes.part, replacement = "")
    genes.part <- gsub(pattern = "var.\\d", x = genes.part, replacement = "")
    return(genes.part)
    }
}


elbow_plot <- function(mat,num_pcs=50,scale=FALSE,center=FALSE,title='',width=3,height=3){
    set.seed(2019)
    mat = data.matrix(mat)
    SVD = irlba(mat, num_pcs, num_pcs,scale=scale,center=center)
    options(repr.plot.width=width, repr.plot.height=height)
    df_plot = data.frame(PC=1:num_pcs, SD=SVD$d);
    p <- ggplot(df_plot, aes(x = PC, y = SD)) +
      geom_point(col="#cd5c5c",size = 1) +
      ggtitle(title)
    return(p)
}

plot.tsne <- function(x, labels,
         main="A tSNE visualization",n=20,
         pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=1) {
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  layout = x

  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=col_vector[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)

  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text,
         col=col_vector[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
return(col_vector[as.integer(labels)])
}

color.vector<-function(x){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    return(col_vector[x])
}

convert.feat.names.to.seurat <- function(feature.names){
    splitted.names <- strsplit(feature.names,split="_")
    return(sapply(splitted.names,function(n){
        paste(gsub(x=n[1],pattern="^chr*",replacement=""),":",n[2],"-",n[3],sep="")
    }))

}

get_plot_limits <- function(plot) {
    gb = ggplot_build(plot)
    xmin = gb$layout$panel_params[[1]]$x.range[1]
    xmax = gb$layout$panel_params[[1]]$x.range[2]
    ymin = gb$layout$panel_params[[1]]$y.range[1]
    ymax = gb$layout$panel_params[[1]]$y.range[2]
    list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

annotate_dotplot_x <- function(g.plot.obj, feat.groups){
    item.2.index <- match(pull(feat.groups,2),unique(pull(feat.groups,2)))
    names(item.2.index) <- pull(feat.groups,1)
    p.limits <- get_plot_limits(g.plot.obj)
    x.items <- as.character(unique(g.plot.obj$data$features.plot))
    x.item.order <- levels(g.plot.obj$data$features.plot)
    x.item.order <- x.item.order[x.item.order %in% x.items]
    n.x.items <- length(x.items)
    x.color.index <- color.vector(item.2.index[x.item.order])
    x.max.pos <- cumsum(rep((p.limits$xmax-p.limits$xmin)/n.x.items,n.x.items))+p.limits$xmin
    x.min.pos <- x.max.pos-(p.limits$xmax-p.limits$xmin)/n.x.items
    color.gene.map<-tibble(colors=x.color.index,geneName=names(item.2.index[x.item.order]))
    data<-cbind(data.frame(xmin=x.min.pos,xmax=x.max.pos,ymin=p.limits$ymin,ymax=p.limits$ymax),join(color.gene.map,feat.groups,by="geneName"))
    gene.type.color.legend<-distinct(data[,c("colors","type")])
    gene.type.color.legend.vector<-gene.type.color.legend[,1]
    names(gene.type.color.legend.vector)<-gene.type.color.legend[,2]
    final.plot <- g.plot.obj + geom_rect(data=data, ymin=0, ymax=p.limits$ymax, aes(xmin=xmin,xmax=xmax,fill=type),colour="white",inherit.aes=F,alpha=0.25) +
    scale_fill_manual(values=gene.type.color.legend.vector,name="Known context")
    return(final.plot)
}


full.snap.to.seurat <- function (obj, eigs.dims = 1:20, norm = TRUE, scale = TRUE)
{
    cat("Epoch: checking input parameters ... \n", file = stderr())
    if (missing(obj)) {
        stop("obj is missing")
    }
    else {
        if (!is.snap(obj)) {
            stop("obj is not a snap object")
        }
        if ((x = nrow(obj)) == 0L) {
            stop("obj is empty")
        }
        if ((x = length(obj@barcode)) == 0L) {
            stop("obj@barcode is empty")
        }
        if ((x = length(obj@file)) == 0L) {
            stop("obj@file is empty")
        }
    }
    if (requireNamespace("Seurat", quietly = TRUE)) {
        require(Seurat)
    }
    else {
        stop("Please install Seurat V3 - learn more at https://github.com/satijalab/seurat")
    }
    if ((x = nrow(obj@gmat)) == 0L) {
        stop("gmat in obj is empty!")
    }
    gmat.use = t(obj@gmat)
    if ((x = nrow(obj@bmat)) > 0L) {
        input.mat = "bmat"
    }
    else if ((x = nrow(obj@pmat)) > 0L) {
        input.mat = "pmat"
    }
    else {
        stop("both pmat and bmat is empty")
    }
    if (input.mat == "bmat") {
        data.use = obj@bmat
        peak.use = as.data.frame(obj@feature)
    }
    else {
        data.use = obj@pmat
        peak.use = as.data.frame(obj@peak)
    }
    if ((x = nrow(data.use)) == 0L) {
        stop("input matrix is empty!")
    }
    metaData.use = obj@metaData
    if ((x = nrow(metaData.use)) == 0L) {
        stop("metaData is empty!")
    }
    ncell = nrow(obj)
    nvar = ncol(obj@smat@dmat)
    if (missing(eigs.dims)) {
        stop("eigs.dims is missing")
    }
    else {
        if (is.null(eigs.dims)) {
            eigs.dims = 1:nvar
        }
        else {
            if (any(eigs.dims > nvar)) {
                stop("'eigs.dims' exceeds PCA dimentions number")
            }
        }
    }
    pca.use = obj@smat@dmat
    if ((x = nrow(pca.use)) == 0L) {
        stop("dimentionality reduction is empty, runLDM first")
    }
    else {
        pca.use = pca.use[, eigs.dims]
    }
    data.use = t(data.use)
    rownames(x = data.use) = peak.use$name
    colnames(x = data.use) = paste0(obj@barcode, 1:ncell)
    colnames(x = gmat.use) = paste0(obj@barcode, 1:ncell)
    rownames(x = pca.use) = paste0(obj@barcode, 1:ncell)
    rownames(metaData.use) = paste0(obj@barcode, 1:ncell)
    pbmc.atac <- CreateSeuratObject(counts = data.use,
        assay = "ATAC")
    pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gmat.use)
    pbmc.atac <- AddMetaData(pbmc.atac, metadata = metaData.use)
    pbmc.atac$tech <- "atac"
    DefaultAssay(pbmc.atac) <- "ATAC"
    colnames(x = pca.use) <- paste0("DC_", eigs.dims)
    pbmc.atac[["SnapATAC"]] <- new(Class = "DimReduc", cell.embeddings = pca.use,
        feature.loadings = matrix(0, 0, 0), feature.loadings.projected = matrix(0,
            0, 0), assay.used = "ATAC", stdev = rep(1, length(eigs.dims)),
        key = "DC_", jackstraw = new(Class = "JackStrawData"),
        misc = list())
    DefaultAssay(pbmc.atac) <- "ACTIVITY"
    if (norm) {
        pbmc.atac <- NormalizeData(pbmc.atac)
    }
    if (scale) {
        pbmc.atac <- ScaleData(pbmc.atac)
    }
    return(pbmc.atac)
}
