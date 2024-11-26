# Some additional functions ----

`%notin%` = Negate(`%in%`)

gene_set_prepare2 <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    paste0(markers_all, collapse=",")
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    paste0(markers_all, collapse=",")
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

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

binarize.expression <- function(seurat.data, assay, grouping=NULL, genes=NULL, cells=NULL, gene.name.correction=NULL, return="bin", cutoff.coef=1.0){
  DefaultAssay(seurat.data) <- assay
  options(dplyr.summarise.inform = FALSE)
  if (!is.null(gene.name.correction)){
    genes[genes %in% names(gene.name.correction)] <- gene.name.correction[names(gene.name.correction)[names(gene.name.correction) %in% genes]] 
  }
  ext.mat <- FetchData(seurat.data, cells = cells, vars = genes)
  if (!is.null(gene.name.correction)){
    colnames(ext.mat)[colnames(ext.mat) %in% gene.name.correction] <- names(gene.name.correction)[gene.name.correction %in% colnames(ext.mat)]
  }
  ident.data <- Idents(seurat.data)
  ext.mat.bin <- apply(ext.mat, 2, function(gene.exp){
    gene.min <- min(gene.exp)
    gene.max <- max(gene.exp)
    if (all(gene.exp==0)){
      gene.exp.bin <- rep(0,length(gene.exp))
      cut.off <- 0
    } else {
      gene.ent <- tibble(gene.exp=gene.exp, ident=ident.data) %>% dplyr::group_by(gr=cut(gene.exp, breaks=seq(gene.min,gene.max*1.1, length.out=25), include.lowest = TRUE), ident) %>%  dplyr::summarise(n=n()) %>% dplyr::group_by(gr) %>% dplyr::mutate(entropy=entropy(n)) %>% dplyr::group_by(gr) %>% dplyr::summarise(entropy=mean(entropy))
      cut.off <- as.numeric(str_extract(as.character(gene.ent[localMinima(gene.ent$entropy)[1],]$gr), pattern="[:digit:]+\\.[:digit:]+"))
      cut.off <- cutoff.coef*cut.off
      gene.exp.bin <- ifelse(gene.exp>cut.off,1,0)
    }
    if (return=="bin"){
      return(gene.exp.bin)
    } else if (return=="cutoff"){
      return(cut.off)
    }
    
  })
  
  if (return=="bin"){
    rownames(ext.mat.bin) <- rownames(ext.mat)
    options(dplyr.summarise.inform = TRUE)
  }
    return(ext.mat.bin)
}

localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

filter.TF <- function(cond.documents, binarized.TF.exp){
  filt.cond.doc <- lapply(names(cond.documents), function(cond.name){
    cluster.number <- str_remove(cond.name, pattern="[:alpha:]{1,2}[:number:]{0,1}\\.")
    if (!is.null(cond.documents[[cond.name]])){
      TF.passed.l.tmp <- lapply(cond.documents[[cond.name]], function(cond){
        TF.expressed.in.condition <- names(which(select(binarized.TF.exp[[as.numeric(cluster.number)]], -starts_with("ident")) %>% colSums() > (nrow(binarized.TF.exp[[as.numeric(cluster.number)]])*.1)))
        TF.passed.l <- cond[str_to_title(str_remove(cond, pattern="_.*")) %in% TF.expressed.in.condition]
        return(TF.passed.l)
      })
    }
    return(TF.passed.l.tmp)
  })
  
  names(filt.cond.doc) <- names(cond.documents)
  return(filt.cond.doc)
}


findpeaks <- function(vec,bw=1,x.coo=c(1:length(vec)))
{
  pos.x.max <- NULL
  pos.y.max <- NULL
  pos.x.min <- NULL
  pos.y.min <- NULL   
  for(i in 1:(length(vec)-1))     {       if((i+1+bw)>length(vec)){
    sup.stop <- length(vec)}else{sup.stop <- i+1+bw
    }
    if((i-bw)<1){inf.stop <- 1}else{inf.stop <- i-bw}
    subset.sup <- vec[(i+1):sup.stop]
    subset.inf <- vec[inf.stop:(i-1)]
    
    is.max   <- sum(subset.inf > vec[i]) == 0
    is.nomin <- sum(subset.sup > vec[i]) == 0
    
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    
    if(is.max & is.nomin){
      pos.x.max <- c(pos.x.max,x.coo[i])
      pos.y.max <- c(pos.y.max,vec[i])
    }
    if(no.max & no.nomin){
      pos.x.min <- c(pos.x.min,x.coo[i])
      pos.y.min <- c(pos.y.min,vec[i])
    }
  }
  return(list(pos.x.max,pos.y.max,pos.x.min,pos.y.min))
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

TF.heatmap <- function(TF.mat.1=NULL,
                       TF.families=NULL,
                       TF.filt=NULL,
                       feature.annotation.region=NULL,
                       cluster.names=NA,
                       links.data=NULL,
                       TF.exprs=FALSE,
                       expr.cutoff=NULL,
                       row.cluster=FALSE,
                       enumerate.features=FALSE,
                       clustering_distance_rows="euclidean",
                       clustering_method_rows="complete") {
  
    TF.used <- rownames(TF.mat.1$per.feat.mat)
    feature.anno <- NULL
    
    if (!is.null(feature.annotation.region)) {
      
      # feature.annotation.gene should be in as per range
      if (is.character(feature.annotation.region)) {
        feature.annotation.region <- StringToGRanges(feature.annotation.region)
      }
      
      # Get all features in region
      features.gr <- StringToGRanges(colnames(TF.mat.1$per.feat.mat))
      
      range_location <- sapply(1:length(features.gr), function (i) {
        feature <- features.gr[i]
        pos <- ifelse(feature %over% feature.annotation.region, "in", "out")
      })
      
      feature.anno <- HeatmapAnnotation(location = as.factor(as.character(range_location)),
                                        col = list(location = c("in" = viridis::cividis(3)[3],
                                                                "out" = viridis::cividis(3)[2])))
    }
    
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
      
      if (enumerate.features) {
        
        number.cols <- sapply(1:ncol(TF.mat.to.plot), function (i) {
          paste0(i, ". ", colnames(TF.mat.to.plot)[i])
        })
        
        colnames(TF.mat.to.plot) <- number.cols
        
      }

      TF.1.plot <- Heatmap(TF.mat.to.plot,
                           cluster_rows = row.cluster,
                           cluster_columns = FALSE,
                           show_row_dend = TRUE,
                           row_names_gp = gpar(fontsize = 6),
                           col = col_fun,
                           row_split = row.split,
                           border = TRUE,
                           row_title_rot = 0,
                           row_gap = unit(2, "mm"),
                           column_names_side = "top",
                           heatmap_legend_param = list(title=cluster.names[1]),
                           bottom_annotation = col_ha,
                           right_annotation = row_ha,
                           top_annotation = feature.anno,
                           clustering_distance_rows = clustering_distance_rows,
                           clustering_method_rows = clustering_method_rows)
      
      return(TF.1.plot)
    } else {
      message("No bindings detected")
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

TF.motifs.per.feature.snakemake <- function(features=NULL, TFBS.data, region, features.in.region=NULL, min.footprint.score=NULL, condition,mc.cores=1){
  # Define features dimension (cols)
  require(parallel)
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

  TFBS.in.features <- mclapply(1:length(TFBS.data), function(x){
    print(paste("Processing",names(TFBS.data)[x],"in",condition,sep=" "))
    tfbs <- TFBS.data[[x]]
    # Subset granges based on bound==1 on given condition
    tfbs.colnames <- colnames(mcols(tfbs))
    tfbs.i <- which(tfbs.colnames==paste(condition,"_bound",sep=""))
    tfbs.filt <- tfbs[elementMetadata(tfbs)[,tfbs.i]==1,]
    tmp.hits <- findOverlaps(query = tfbs.filt, subject = features.in.region, minoverlap = 1)
    tmp.features <- GRangesToString(features.in.region[subjectHits(tmp.hits)])
    tfbs.hits <- tfbs.filt[queryHits(tmp.hits)]
    if (length(tfbs.hits)>0){
      tfbs.hits$features.with.hits <- tmp.features
    }
    return(tfbs.hits)
  }, mc.cores=mc.cores)

  names(TFBS.in.features) <- names(TFBS.data)
  TF.hit.count <- sapply(TFBS.in.features, length)
  TF.hits <- TFBS.in.features[TF.hit.count>0]

  print(paste("Found total of ", sum(TF.hit.count>0), " TF hits the region",sep=""))

  # Create zero matrix
  TF.motif.matrix <- matrix(0, nrow = length(TFBS.data), ncol=length(features.in.region))
  rownames(TF.motif.matrix) <- names(TFBS.data)
  colnames(TF.motif.matrix) <- GRangesToString(features.in.region)

  TF.hit.coordinates <-  lapply(TFBS.in.features, function(x){
      tmp <- start(x)
      names(tmp) <- x$features.with.hits
      return(tmp)
    })
  
  names(TF.hit.coordinates) <- str_extract(names(TFBS.in.features) ,pattern = ".*[[:digit:]]{1}\\.[[:alpha:]]{1}_") %>% str_sub(start=1, end=-2)
    
  # Loop over all TFBS binding events which overlapped features in the gene region
  lapply(names(TF.hits), function(tf){
    features.with.hits <- TF.hits[[tf]]$features.with.hits
    tfbs.colnames <- colnames(mcols(TF.hits[[tf]]))
    tfbs.i <- which(tfbs.colnames==paste(condition,"_score",sep=""))
    footprint.scores <- TF.hits[[tf]][,tfbs.i]

    avg.footprint.score.per.feat <- tapply(INDEX=features.with.hits, X=mcols(footprint.scores)[,1], FUN=mean)
    TF.motif.matrix[tf,names(avg.footprint.score.per.feat)] <<- avg.footprint.score.per.feat
  })
  return(list(per.feat.mat=TF.motif.matrix, TF.hit.coordinates=TF.hit.coordinates))
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

TF.motifs.per.feature.snakemake_v2 <- function(features=NULL, TFBS.data, region, features.in.region=NULL, min.footprint.score=NULL, condition,mc.cores=1){
  # Define features dimension (cols)
  require(parallel)
  require(valr)
  if (!is.null(features)){
    features.gr <- StringToGRanges(features)
    names(features.gr) <- rownames(features)
    features.in.region <- features.gr[features.gr %over% StringToGRanges(region)]
  } else {
    features.in.region <- StringToGRanges(features.in.region)
  }
  #print(paste("Found ", length(features.in.region), " features in the region", sep=""))
  if (!class(TFBS.data)=="CompressedGRangesList"){
    TFBS.data <- GRangesList(TFBS.data)
  }
  
  TFBS.in.features <- mclapply(1:length(TFBS.data), function(x){
    print(paste("Processing",names(TFBS.data)[x],"in",condition,sep=" "))
    tfbs <- TFBS.data[[x]]
    # Subset granges based on bound==1 on given condition
    tfbs.colnames <- colnames(mcols(tfbs))
    tfbs.i <- which(tfbs.colnames==paste(condition,"_bound",sep=""))
    tfbs.filt <- tfbs[elementMetadata(tfbs)[,tfbs.i]==1,]
    tmp.hits <- findOverlaps(query = tfbs.filt, subject = features.in.region, minoverlap = 1)
    tmp.features <- GRangesToString(features.in.region[subjectHits(tmp.hits)])
    tfbs.hits <- tfbs.filt[queryHits(tmp.hits)]
    if (length(tfbs.hits)>0){
      tfbs.hits$features.with.hits <- tmp.features
    }
    return(tfbs.hits)
  }, mc.cores=mc.cores)
  
  names(TFBS.in.features) <- names(TFBS.data)
  TF.hits <- sapply(1:length(TFBS.in.features), function(x){
    tmp <- gr_to_bed(TFBS.in.features[[x]])
    if (nrow(tmp)>0){
      tibble(tmp,motif_id=names(TFBS.in.features)[x])
    }
  })
  
  return(bind_rows(TF.hits))  

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

precalculate.matches <- function(db.conn, table.name, features.gr){
  require(magrittr)
  require(RSQLite)
  require(valr)
  require(GenomicRanges)
  conn <- dbConnect(RSQLite::SQLite(), db.conn)
  table.name <- table.name
  features.tb <- gr_to_bed(features.gr)
  features.strings <- paste(pull(features.tb, chrom),"-",pull(features.tb, start),"-",pull(features.tb, end),sep="")
  
  snakemake <- dbGetQuery(conn, paste('SELECT TFBS_chr, TFBS_start, TFBS_end, id FROM ',mydb.table, ';',sep=""))
  colnames(snakemake) <- c("chr","start","end","id")
  snakemake.gr <- GenomicRanges::makeGRangesFromDataFrame(snakemake)
  snake.hits <- findOverlaps(query = snakemake.gr, subject = features.gr, minoverlap = 1L)
  
  snake.hits.tb <- tibble(motif_id=pull(snakemake, id)[queryHits(snake.hits)], feature=features.strings[subjectHits(snake.hits)])
  dbWriteTable(conn, name=table.name, value=snake.hits.tb, append = TRUE) 
}


RenameGenesSeurat <- function(obj, newnames,assay="RNA") { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  obj[['RNA_name']] <- obj[[assay]]
  RNA <- obj@assays$RNA_name
  if (length(RNA@scale.data) > 0){
    tmp.conv <- tibble(id=RNA@counts@Dimnames[[1]], symbol=newnames)
  }

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts) >0 & class(RNA@data)[1]=="dgCMatrix") {RNA@counts@Dimnames[[1]]            <- newnames}
    if (length(RNA@data) >0 ){ RNA@data@Dimnames[[1]]                <- newnames}
    if (length(RNA@scale.data) > 0 & !is.matrix(RNA@scale.data)){RNA@scale.data@Dimnames[[1]]    <- newnames}
    if (length(RNA@scale.data) > 0 & is.matrix(RNA@scale.data)){rownames(RNA@scale.data)    <- tmp.conv$symbol[match(rownames(RNA@scale.data),tmp.conv$id)]}
    #if (length(RNA@scale.data)) dimnames(RNA@scale.data)[[1]]    <- tmp.conv$symbol[match(dimnames(RNA@scale.data)[[1]],tmp.conv$id)]
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA_name <- RNA
  return(obj)
}


RenameGenesSeurat_v2 <- function(obj, newnames,assay="RNA") { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  obj[['RNA_name']] <- obj[[assay]]
  RNA <- obj@assays$RNA_name
  if (length(RNA@layers[["scale.data"]]) > 0){
    tmp.conv <- tibble(id=RNA@layers[["counts"]]@Dimnames[[1]], symbol=newnames)
  }
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@layers[["counts"]]) >0 & class(RNA@layers[["data"]])[1]=="dgCMatrix") {RNA@layers[["counts"]]@Dimnames[[1]]<- newnames}
    if (length(RNA@layers[["data"]]) >0 ){RNA@layers[["data"]]@Dimnames[[1]] <- newnames}
    if (length(RNA@layers[["scale.data"]]) > 0 & !is.matrix(RNA@layers[["scale.data"]])){RNA@layers[["scale.data"]]@Dimnames[[1]]    <- newnames}
    if (length(RNA@layers[["scale.data"]]) > 0 & is.matrix(RNA@layers[["scale.data"]])){rownames(RNA@layers[["scale.data"]])    <- tmp.conv$symbol[match(rownames(RNA@layers[["scale.data"]]),tmp.conv$id)]}
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
    message("Warning: Found non-ENS ID for argument feature format 'ens'")
  }
  if (feature.format == "symbol" & any(grepl("*ENS", features))) {
    message("Warning: Found ENS ID for argument feature format 'symbol'")
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

    succ.matches <- sum(!is.na(v.out))
    n.tot <- length(features)
    
    msg <- str_interp("Instance: Found matching for ${succ.matches} features out of total ${n.tot} provided features")
    
    message(msg)

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


ChromVar_DA <- function(da.features, id.1, id.2, motifs, seurat.object, genome = BSgenome.Mmusculus.UCSC.mm10){
  
  da_features[[id.1]][[id.2]] %>% dplyr::filter((avg_log2FC > 0.5 | avg_log2FC < -0.5) & p_val_adj < 0.05) %>% pull(feature) -> sub.da.features
  
  tmp.1.motifs <- CreateMotifMatrix(
    features = StringToGRanges(sub.da.features),
    pwm = motifs,
    genome = genome,
    score = FALSE,
    use.counts = FALSE,
    sep = c("-", "-")
  )
  
  seurat.object <- RunChromVAR(
    object = subset(seurat.object, features=sub.da.features),
    motif.matrix = tmp.1.motifs,
    assay="peaks",
    genome = BSgenome.Mmusculus.UCSC.mm10,
    new.assay.name = "chromvar.da"
  )
  
  DefaultAssay(seurat.object) <- "chromvar.da"
  markers_chromvar <- FindMarkers(
    object = seurat.object,
    only.pos = FALSE,
    ident.1 = id.1,
    ident.2 = id.2,
    test.use = 'LR',
    latent.vars = 'nCount_peaks'
  ) %>% dplyr::filter(p_val_adj <= 0.05 & (avg_log2FC >= 0.5 | avg_log2FC <= -0.5)) %>% rownames_to_column(var="motif.id") %>% as_tibble() %>% arrange(avg_log2FC)
  
  return(markers_chromvar)

# Complement of %in% operator
`%notin%` <- purrr::negate(`%in%`)


# Version from 12.08.2022
group.TF.feature.heatmap <- function(mat.list, ident.names, feature.of.interest){
  
  #' --------------------------------------------------------------------------
  #' @param mat.list  A named list of matrices
  #'---------------------------------------------------------------------------
  
  # NOTE: Set3 contains 12 colors.
  # RColorBrewer::brewer.pal rotates the set if index if larger
  n.matrices <- length(mat.list)
  
  # Compute union over row names of the matrix list
  motif.union <- lapply(mat.list, function (mat) rownames(mat)) %>% unlist() %>% unique()
  
  # Transform each matrix into same motif space
  mat.list.full <- lapply(mat.list, function (mat) {
    
    exterior.motifs <- motif.union[motif.union %notin% rownames(mat)]
    
    zero.mat <- matrix(0, length(exterior.motifs), ncol(mat))
    
    rownames(zero.mat) <- exterior.motifs
    colnames(zero.mat) <- colnames(mat)
    
    out.mat <- rbind(mat, zero.mat)
    
    rownames(out.mat) <- str_replace(rownames(out.mat), "(.*)\\.[A-Z]_", "")
    
    out.mat <- out.mat[order(row.names(out.mat)), ]
    
    out.f <- out.mat[,feature.of.interest]
    
    return (out.f)
  })
  
  joined.feature.tf.mat <- do.call(cbind, mat.list.full)
  
  colnames(joined.feature.tf.mat) <- ident.names
  
  hclust.out <- dist(joined.feature.tf.mat) %>% hclust()
  
  joined.feature.tf.mat <- joined.feature.tf.mat[hclust.out$order,]
  
  joined.feature.tf.mat <- joined.feature.tf.mat[apply(joined.feature.tf.mat[,-1], 1, function(x) !all(x == 0)), ]
  
  hm.out <- Heatmap(joined.feature.tf.mat, cluster_rows = F, cluster_columns = F,
                    name = feature.of.interest, col = viridis::viridis(10))
  return(hm.out)
}}

read_consensus_bed <- function(bed_file){
  data.tmp <- read.table(bed_file,sep="\t",fill=TRUE,header = F,col.names=c("chr", "start", "end", "replicate_peak_start_coordinates","replicate_peak_end_coordinates","total_signal","max_signal","max_signal_region","replicate","consensus_count"))
  return(data.tmp)
}

overlap.bounds <- function (range, fpath, res.type = "bound", col.names = "Default", TF = NULL)
{
  #'
  #' @param range     (str)       Range over which to visualize results
  #' @param fpath     (str)       Path to TOBIAS BINDetect results
  #' @param res.type  (str)       Type of such results "bound"/"unbound"
  #' @param col.names (str)       Column names, default "Default
  #' @param TF        (str)       TF to which restrict results
  #'
  
  # Dependencies
  require("Signac")
  require("GenomicRanges")
  
  if (is.character(range)) {
    range <- StringToGRanges(range)
  }
  
  bound <- ifelse(res.type == "bound", T, F)
  
  bd.res <- get_BINDetect_results(fpath, bound = bound, col.names = col.names)
  
  overlaps <- lapply(bd.res, function (r) {
    out <- r[r %over% range]
  })
  
  if (!is.null(TF)) {
    overlaps <- overlaps[names(overlaps) == TF][[1]]
    overlaps <- overlaps[,c(1,2,3)]
  }
  
  return (overlaps)
  
}

consensus_bed2Granges <- function(bed_file){
  require(GenomicRanges)
  data_tmp <- read_consensus_bed(bed_file)
  consensus.bed.granges <- makeGRangesFromDataFrame(data_tmp,
                                                   keep.extra.columns = TRUE,
                                                   seqnames.field = "chr",
                                                   start.field = "start",
                                                   end.field = "end",
                                                   ignore.strand = T)
  return(consensus.bed.granges)
}



findEChO <- function(EChO.matrix, span, foci, EChO.thr=120){
  coordinates <- apply(foci[span$line,] %>% select(c("chr","start","end")), 1, function(r){r<-str_remove_all(string = r,pattern = " ");paste(r,collapse = "-")})
  rownames(EChO.matrix) <- coordinates
  #makeGRangesFromDataFrame(foci[span$line,] %>% select(c("chr","start","end")))
  EChO.true.i <- which(rowSums(EChO.matrix<120)>0)
  return(list(EChO.true.i=EChO.true.i,coordinates=coordinates[EChO.true.i]))
}

plotSmoothedAccessibility_v2 <- function(dataset, covariate.genes.to.plot, features.df, gene.feat.name,gene_id2name, gene_name2id, rV2=TRUE, drop.non.linked=FALSE){
  
  # Map gene symbols to Ensembl IDs using the hash
  covariate.ids.to.plot <- sapply(covariate.genes.to.plot, function(f) { gene_name2id[[f]] })
  
  # Check for NAs (genes not found in the mapping)
  if (any(is.na(covariate.ids.to.plot))) {
    missing_genes <- covariate.genes.to.plot[is.na(covariate.ids.to.plot)]
    stop("The following genes were not found in gene_name2id mapping: ", paste(missing_genes, collapse = ", "))
  }
  
  # Ensure covariate.ids.to.plot is a character vector
  covariate.ids.to.plot <- as.character(covariate.ids.to.plot)
  
  # Fetch peaks data
  DefaultAssay(dataset) <- "peaks"
  features.to.show <- features.df$feature
  peaks.data <- FetchData(dataset, vars = features.to.show)
  features.granges.gr <- StringToGRanges(features.to.show)
  
  # Fetch expression data using Ensembl IDs
  DefaultAssay(dataset) <- "RNA"
  expression.data <- FetchData(dataset, vars = covariate.ids.to.plot)
  
  # Assign gene symbols as column names
  colnames(expression.data) <- covariate.genes.to.plot
  
  # Prepare covariate data
  eoi <- data.frame(expression.data)
  eoi$pseudotime <- dataset$VIA_pseudotime
  eoi <- eoi %>% arrange(pseudotime)
  
  # Prepare peaks data
  peaks.data.ss <- peaks.data
  peaks.data.ss$pseudotime <- dataset$VIA_pseudotime
  peaks.data.ss$label <- dataset$rv2.lineage
  peaks.data.ss$barcode <- rownames(peaks.data.ss)
  peaks.data.ss <- peaks.data.ss %>% arrange(pseudotime)
  
  # Add covariate gene expression to peaks data
  peaks.data.ss <- cbind(peaks.data.ss, eoi[covariate.genes.to.plot])
  
  # Calculate rolling mean for covariate genes
  covariate.rolled.means <- peaks.data.ss %>%
    dplyr::select(all_of(covariate.genes.to.plot)) %>%
    zoo::rollapply(width = 6, by = 1, FUN = mean, align = "center", by.column = TRUE, fill = c(NA, NA, NA))
  colnames(covariate.rolled.means) <- paste(colnames(covariate.rolled.means), "_rolled", sep = "")
  peaks.data.ss <- cbind(peaks.data.ss, covariate.rolled.means)
  
  # Define color mapping for pseudotime
  cols_pseudotime <- list(pseudotime = viridis::viridis(length(unique(peaks.data.ss$pseudotime))))
  names(cols_pseudotime$pseudotime) <- unique(peaks.data.ss$pseudotime)
  
  # Prepare row labels with adjusted TSS distance
  feature_labels <- sapply(features.df$Gene_TSS_distance, function(dist){
    if (is.na(dist)) {
      return(NA)
    } else if (dist < 0) {
      paste0(gene.feat.name," +", round(abs(dist) / 1000, 1), "kb")
    } else {
      paste0(gene.feat.name," -", round(abs(dist) / 1000, 1), "kb")
    }
  })
  
  # Ensure the feature labels correspond to the order of features.to.show
  names(feature_labels) <- features.df$feature
  feature_labels <- feature_labels[features.to.show]
  
  # Adjust data based on rV2 value
  if (rV2){
    # rV2 is TRUE
    # Include PRO1, PRO2, CO1, CO2, and all GA groups
    ga_labels <- grep("^GA", unique(peaks.data.ss$label), value = TRUE)
    included_labels <- c("PRO1", "PRO2", "CO1", "CO2", ga_labels)
  } else {
    # rV2 is FALSE
    # Include PRO1, PRO2, CO1, CO2, and all GL groups
    gl_labels <- grep("^GL", unique(peaks.data.ss$label), value = TRUE)
    included_labels <- c("PRO1", "PRO2", "CO1", "CO2", gl_labels)
  }
  
  # Subset the data for included labels
  data_subset <- peaks.data.ss %>% filter(label %in% included_labels)
  
  # Ensure labels are factors with levels in the desired order
  data_subset$label <- factor(data_subset$label, levels = included_labels)
  
  # Scale and smooth data
  data.scaled <- dplyr::select(data_subset, all_of(features.to.show)) %>%
    as.matrix() %>% scale(scale = TRUE, center = TRUE)
  data.scaled[!is.finite(data.scaled)] <- 0
  data.scaled <- apply(data.scaled, MARGIN = 2, FUN = function(vec){
    smoo <- smooth.spline(vec, cv = FALSE, penalty = 0.8)
    smoo$y
  })
  rownames(data.scaled) <- data_subset$barcode
  
  # Prepare annotations
  labels <- data_subset$label
  
  ha.bot <- HeatmapAnnotation(
    pseudotime = data_subset$pseudotime,
    col = cols_pseudotime,
    show_legend = FALSE,
    annotation_name_gp = grid::gpar(fontsize = 20),
    simple_anno_size = unit(2, "cm")
  )
  
  cov.data <- dplyr::select(data_subset, ends_with("_rolled"))
  colnames(cov.data) <- covariate.genes.to.plot
  expression_col_fun <- lapply(seq_along(covariate.genes.to.plot), function(i){
    gene <- covariate.genes.to.plot[i]
    colorRamp2(c(min(cov.data[[i]], na.rm = TRUE), max(cov.data[[i]], na.rm = TRUE)), c("white", "orange"))
  })
  names(expression_col_fun) <- covariate.genes.to.plot
  
  ha.top <- HeatmapAnnotation(
    df = cov.data,
    annotation_label = covariate.genes.to.plot,
    col = expression_col_fun,
    simple_anno_size = unit(3, "cm"),
    height = unit(8, "cm"),
    annotation_name_gp = grid::gpar(fontsize = 20),
    gp = grid::gpar(fontsize = 20)
  )
  
  # Create heatmap
  p <- Heatmap(
    t(data.scaled) + 1,
    name = "Accessibility",
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_labels = feature_labels,
    bottom_annotation = ha.bot,
    top_annotation = ha.top,
    cluster_rows = FALSE,
    column_split = labels,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    col = viridis::magma(100),
    use_raster = TRUE,
    column_names_gp = grid::gpar(fontsize = 20),
    row_names_gp = grid::gpar(fontsize = 15),
    column_title_gp = grid::gpar(fontsize = 20)
  )
  
  return(p)
}
  
plotSmoothedAccessibility <- function(dataset, covariate.genes.to.plot, region.of.interest, rV2=TRUE, drop.non.linked=FALSE){
  
  DefaultAssay(dataset) <- "peaks"
  peaks.data <- FetchData(dataset, vars = rownames(dataset))
  all.features <- rownames(dataset)
  all.features.gr <- StringToGRanges(all.features)
  DefaultAssay(dataset) <- "RNA"
  expression.data <- FetchData(dataset, vars = rownames(dataset))
  
  features.granges.gr <- all.features.gr[all.features.gr %over% StringToGRanges(region.of.interest)]
  features.in.region <- GRangesToString(features.granges.gr)
  
  DefaultAssay(dataset) <- "peaks"
  f1.range <- range(Annotation(dataset)[Annotation(dataset)$gene_name %in% covariate.genes.to.plot[1]], ignore.strand = TRUE)
  f2.range <- f1.range+50000
  
  # dplyr::select features to be included: coding+promixal+significant distal ones
  links.all <- Links(dataset)[Links(dataset) %over% features.granges.gr & elementMetadata(Links(dataset))$pvalue < 0.01]
  links.all.peaks.gr <- StringToGRanges(elementMetadata(links.all)$peak)
  
  features.to.show.gr <- sort(unique(c(features.granges.gr[features.granges.gr %over% f2.range],links.all.peaks.gr)))
  features.to.show <- GRangesToString(features.to.show.gr)
  
  links.to.highlight <- Links(dataset)[Links(dataset) %over% f2.range & elementMetadata(Links(dataset))$pvalue < 0.01]
  links.to.highlight.peak.gr <- StringToGRanges(elementMetadata(links.to.highlight)$peak)
  
  if (drop.non.linked){
    features.to.show.gr <- features.to.show.gr[features.to.show.gr %over% links.to.highlight.peak.gr]
    features.to.show <- GRangesToString(features.to.show.gr)
  }
 
  # Preparing covariate data
  id.to.plot <- convert_feature_identity(dataset, "RNA", covariate.genes.to.plot, "symbol")
  eoi <- data.frame(expression.data[,colnames(expression.data) %in% id.to.plot])
  colnames(eoi) <- convert_feature_identity(dataset, "RNA", colnames(eoi), "ens")
  rownames(eoi) <- rownames(expression.data)
  eoi$pseudotime <- dataset$VIA_pseudotime
  eoi <- eoi %>% arrange(pseudotime)

  peaks.data.ss <- peaks.data[,features.to.show]
  peaks.data.ss <- merge(peaks.data.ss, data.frame(dataset$VIA_pseudotime), by = "row.names")
  rownames(peaks.data.ss) <- peaks.data.ss$Row.names
  peaks.data.ss <- merge(peaks.data.ss, data.frame(dataset$rv2.lineage), by = "row.names") 
  rownames(peaks.data.ss) <- peaks.data.ss$Row.names
  peaks.data.ss <- peaks.data.ss[,3:ncol(peaks.data.ss)]
  col.names <- colnames(peaks.data.ss)
  col.names[length(col.names) - 1] <- "pseudotime"
  col.names[length(col.names)] <- "label"
  colnames(peaks.data.ss) <- col.names
 
  f1.expression <- data.frame(eoi %>% dplyr::select(all_of(covariate.genes.to.plot)), row.names = rownames(eoi))
  peaks.data.ss <- merge(peaks.data.ss, f1.expression, by = "row.names")
  rownames(peaks.data.ss) <- peaks.data.ss$Row.names
  peaks.data.ss <- peaks.data.ss[,2:ncol(peaks.data.ss)] 
  
  peaks.data.ss <- peaks.data.ss %>% arrange(pseudotime, id.to.plot[1])
  
  covariate.rolled.means <- rownames_to_column(peaks.data.ss, var = "barcode") %>% as_tibble() %>% dplyr::select(contains(str_replace(string = covariate.genes.to.plot, pattern = "-", replacement = "."))) %>% zoo::rollapply(width = 6, by = 1, FUN = mean, align = "center", by.column=TRUE, fill=c(0))
  colnames(covariate.rolled.means) <- paste(colnames(covariate.rolled.means),"_rolled",sep="")
  peaks.data.ss <- cbind(peaks.data.ss,covariate.rolled.means)
  
  cols_pseudotime <- list(pseudotime = viridis::viridis(length(unique(peaks.data.ss$pseudotime))))
  names(cols_pseudotime$pseudotime) <- unique(peaks.data.ss$pseudotime)

  data.scaled.tmp <- dplyr::select(peaks.data.ss, starts_with("chr")) %>% as.matrix() %>% scale(scale = T, center = T) 
  feature.filter.l <- colSums(is.finite(data.scaled.tmp))>0
  data.scaled.tmp[!is.finite(data.scaled.tmp)] <- 0
  data.scaled <-apply(data.scaled.tmp,MARGIN=2,FUN = function(vec){smoo <- smooth.spline(vec, cv = F, penalty = 0.8);smoo$y})
  rownames(data.scaled) <- rownames(peaks.data.ss)
  
  # if (any(is.nan(data.scaled))) {
  #   data.scaled[is.nan(data.scaled)] <- 0
  # }

  # Defining rowannotation colors
  # DefaultAssay(dataset) <- "peaks"
  # features.granges <- StringToGRanges(features.to.show)
  # f1.range <- range(Annotation(dataset)[Annotation(dataset)$gene_name %in% covariate.genes.to.plot[1]], ignore.strand = TRUE)
  # f2.range <- f1.range+50000
  #features.in.f1.range <- GRangesToString(features.granges[features.granges %over% f1.range])
  
  range_location <- sapply(1:length(features.to.show.gr), function (i) {
    feature <- features.to.show.gr[i]
    case_when(feature %over% f1.range ~ "coding",
                     feature %over% f2.range & !(feature %over% f1.range) ~ "proximal",
                     !(feature %over% f2.range) & !(feature %over% f1.range) ~ "distal")
  })
  
  ha_row_left <- rowAnnotation(
    location = as.factor(as.character(range_location)), 
    col = list(location = c("coding" = "firebrick2","proximal"="cornflowerblue","distal"="lightskyblue")),
    annotation_name_gp = grid::gpar(fontsize = 20),
    simple_anno_size = unit(2, "cm")
  )
  
  mark_at <- which(features.to.show.gr %over% links.to.highlight.peak.gr)
  
  ha_row_right <- rowAnnotation(
    Marked_features = anno_mark(at=mark_at,labels=GRangesToString(links.to.highlight.peak.gr)), 
    annotation_name_gp = grid::gpar(fontsize = 20),
    simple_anno_size = unit(2, "cm")
  )
  
  # k <- ceiling(9 / length(features.in.region)) # 9 == number of unique colors
  # cols.named <- rep(RColorBrewer::brewer.pal(9, "Set1"), k)[1:length(features.in.region)]
  # names(cols.named) <- features.in.region
  # 
  # ha_right <- rowAnnotation(
  #   peak = features.in.region,
  #   col = list(peak = cols.named),
  #   show_legend = F
  # )
  
  # Preparing data matrix itself
  #heatmap.data <- as.matrix(peaks.data.ss[1:(ncol(peaks.data.ss) - (2 + 2 * length(covariate.genes.to.plot)))]) # (pseudotime + label) + (gene + gene.rolled) * n(genes)
  fac.ord <- aggregate(peaks.data.ss$pseudotime, list(peaks.data.ss$label), mean) %>% arrange(x) %>% pull(Group.1)
  
  row.anno <- sapply(colnames(data.scaled), function (id) {
    f.start <- strsplit(id, "-")[[1]][2]
    f.end <- strsplit(id, "-")[[1]][3]
    
    as.numeric(f.end) - as.numeric(f.start)
  })
  
  if (rV2){
  cells.non.gl <- Cells(subset(dataset, subset = rv2.lineage %notin% sapply(1:5, function(n) {paste0("GL", n)})))
  data.scaled.ga <- data.frame(data.scaled) %>% filter(rownames(data.scaled) %in% cells.non.gl)
  colnames(data.scaled.ga) <- str_replace_all(colnames(data.scaled.ga), "\\.", "-")
  
  cells.non.ga <- Cells(subset(dataset, subset = rv2.lineage %notin% sapply(1:6, function(n) {paste0("GA", n)})))
  data.scaled.gl <- data.scaled[rownames(data.scaled) %in% cells.non.ga,]
  colnames(data.scaled.gl) <- str_replace_all(colnames(data.scaled.gl), "\\.", "-")
  
  ha.bot.ga <- HeatmapAnnotation(pseudotime = peaks.data.ss[rownames(peaks.data.ss) %in% cells.non.gl,]$pseudotime, 
                                 col = cols_pseudotime, show_legend = F, annotation_name_gp = grid::gpar(fontsize = 20),
                                 simple_anno_size = unit(2, "cm"))
  labels.ga <- factor(peaks.data.ss[rownames(peaks.data.ss) %in% cells.non.gl,]$label, levels = fac.ord)
  
  ha.bot.gl <- HeatmapAnnotation(pseudotime = peaks.data.ss[rownames(peaks.data.ss) %in% cells.non.ga,]$pseudotime, 
                                 col = cols_pseudotime, show_legend = F, annotation_name_gp = grid::gpar(fontsize = 20),
                                 simple_anno_size = unit(2, "cm"))
  labels.gl <- factor(peaks.data.ss[rownames(peaks.data.ss) %in% cells.non.ga,]$label, levels = fac.ord)

  cov.gl.data <- dplyr::select(peaks.data.ss[rownames(peaks.data.ss) %in% cells.non.ga,], contains("rolled"))
  colnames(cov.gl.data) <- covariate.genes.to.plot
  gl.expression_col_fun = replicate(expr = colorRamp2(c(min(cov.gl.data), max(cov.gl.data)), c("white", "orange")),n=ncol(cov.gl.data))
  names(gl.expression_col_fun) <- covariate.genes.to.plot
  
  cov.ga.data <- dplyr::select(peaks.data.ss[rownames(peaks.data.ss) %in% cells.non.gl,], contains("rolled"))
  colnames(cov.ga.data) <- covariate.genes.to.plot
  ga.expression_col_fun = replicate(expr = colorRamp2(c(min(cov.ga.data), max(cov.ga.data)), c("white", "orange")),n=ncol(cov.ga.data))
  names(ga.expression_col_fun) <- covariate.genes.to.plot
  
  ha.top.ga <- HeatmapAnnotation(df=cov.ga.data,
                                 annotation_label = covariate.genes.to.plot,
                                 col= ga.expression_col_fun ,
                                 simple_anno_size = unit(3, "cm"), height = unit(8, "cm"),
                                 annotation_name_gp = grid::gpar(fontsize = 20),
                                 gp = grid::gpar(fontsize = 20))
  
  ha.top.gl <- HeatmapAnnotation(df=cov.gl.data,
                                 annotation_label = covariate.genes.to.plot,
                                 col= gl.expression_col_fun ,
                                 simple_anno_size = unit(3, "cm"), height = unit(8, "cm"),
                                 annotation_name_gp = grid::gpar(fontsize = 20),
                                 gp = grid::gpar(fontsize = 20))
  
  col.p1 <- colorRamp2(breaks=c(seq(min(data.scaled.ga)+1,max(data.scaled.gl)+1,length.out=100)),colors=viridis::magma(100))

  p.1 <- Heatmap(t(data.scaled.ga) + 1,
                 show_column_names = F, 
                 show_row_names = F,
                 bottom_annotation = ha.bot.ga,
                 top_annotation = ha.top.ga,
                 left_annotation = ha_row_left,
                 right_annotation = ha_row_right,
                 cluster_rows = F,
                 column_split = labels.ga,
                 cluster_columns = F,
                 cluster_column_slices = F,
                 #col = viridis::magma(100),
                 col=col.p1,
                 row_labels = paste(1:length(colnames(data.scaled.ga)), ".\nLength: ", row.anno, sep = ''),
                 use_raster = T,
                 column_names_gp = grid::gpar(fontsize = 20),
                 row_names_gp = grid::gpar(fontsize = 20),
                 column_title_gp = grid::gpar(fontsize = 20))
  
  col.p2 <- colorRamp2(breaks=c(seq(min(data.scaled.gl)+1,max(data.scaled.gl)+1,length.out=100)),colors=viridis::magma(100))
  
  p.2 <- Heatmap(t(data.scaled.gl) + 1,
                 show_column_names = F,
                 show_row_names = F,
                 bottom_annotation = ha.bot.gl,
                 top_annotation = ha.top.gl,
                 left_annotation = ha_row_left,
                 right_annotation = ha_row_right,
                 cluster_rows = F,
                 column_split = labels.gl,
                 cluster_columns = F,
                 cluster_column_slices = F,
                 col = col.p2,
                 row_labels = paste(1:length(colnames(data.scaled.gl)), ".\nLength: ", row.anno, sep = ''),
                 use_raster = T,
                 column_names_gp = grid::gpar(fontsize = 20),
                 row_names_gp = grid::gpar(fontsize = 20),
                 column_title_gp = grid::gpar(fontsize = 20))

  return(list(p.1,p.2))
  } else {
    colnames(data.scaled) <- str_replace_all(colnames(data.scaled), "\\.", "-")

    ha.bot <- HeatmapAnnotation(pseudotime = peaks.data.ss$pseudotime, 
                                   col = cols_pseudotime, show_legend = F, annotation_name_gp = grid::gpar(fontsize = 20),
                                   simple_anno_size = unit(2, "cm"))
    labels <- factor(peaks.data.ss$label, levels = fac.ord)
    
    ha.bot <- HeatmapAnnotation(pseudotime = peaks.data.ss$pseudotime, 
                                   col = cols_pseudotime, show_legend = F, annotation_name_gp = grid::gpar(fontsize = 20),
                                   simple_anno_size = unit(2, "cm"))
    labels <- factor(peaks.data.ss$label, levels = fac.ord)
    
    cov.data <- dplyr::select(peaks.data.ss, contains("rolled"))
    colnames(cov.data) <- covariate.genes.to.plot
    expression_col_fun = replicate(expr = colorRamp2(c(min(cov.data), max(cov.data)), c("white", "orange")),n=ncol(cov.data))
    names(expression_col_fun) <- covariate.genes.to.plot
    
    ha.top <- HeatmapAnnotation(df=cov.data,
                                   annotation_label = covariate.genes.to.plot,
                                   col= expression_col_fun ,
                                   simple_anno_size = unit(3, "cm"), height = unit(8, "cm"),
                                   annotation_name_gp = grid::gpar(fontsize = 20),
                                   gp = grid::gpar(fontsize = 20))
    

    p.1 <- Heatmap(t(data.scaled) + 1,
                   show_column_names = F, 
                   show_row_names = F,
                   bottom_annotation = ha.bot,
                   top_annotation = ha.top,
                   left_annotation = ha_row_left,
                   right_annotation = ha_row_right,
                   cluster_rows = F,
                   column_split = labels,
                   cluster_columns = F,
                   cluster_column_slices = F,
                   col = viridis::magma(100),
                   row_labels = paste(1:length(colnames(data.scaled)), ".\nLength: ", row.anno, sep = ''),
                   use_raster = T,
                   column_names_gp = grid::gpar(fontsize = 20),
                   row_names_gp = grid::gpar(fontsize = 20),
                   column_title_gp = grid::gpar(fontsize = 20))
    return(p.1)
  }
}

formatGOenrichRes <- function(sampleGOdata,p.val.thr=0.05){
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  allRes <- as_tibble(GenTable(sampleGOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 5000))
  specRes <- allRes %>% filter(GO.ID %in% c("GO:0097154","GO:1905962","GO:0045664"))
  allRes <- allRes %>% filter(classicFisher < p.val.thr)
  
  allRes$genes.in.classes <- sapply(allRes$GO.ID,function(GO_id){
    tmp.genes <- intersect(genesInTerm(sampleGOdata, GO_id)[[1]], sigGenes(sampleGOdata))
    return(paste(tmp.genes, collapse = ", "))
  })
  
  specRes$genes.in.classes <- sapply(specRes$GO.ID,function(GO_id){
    tmp.genes <- intersect(genesInTerm(sampleGOdata, GO_id)[[1]], sigGenes(sampleGOdata))
    return(paste(tmp.genes, collapse = ", "))
  })
  
  return(list(allRes=allRes,specRes=specRes))
}


ExpandAxes <- function(x, y, expansion_factor = 0.1) {
  x_range <- range(x)
  y_range <- range(y)
  
  # Expand the range by a given factor (default is 10%)
  x_expand <- diff(x_range) * expansion_factor
  y_expand <- diff(y_range) * expansion_factor
  
  return(list(x.range = c(x_range[1] - x_expand, x_range[2] + x_expand),
              y.range = c(y_range[1] - y_expand, y_range[2] + y_expand)))
}

library(dplyr)
library(ComplexHeatmap)
library(circlize)  # Needed for color scaling in ComplexHeatmap

# Define the function to compute hypergeometric p-values between scRNA and scATAC clusters and return heatmap
hypergeometric_test_clusters <- function(scRNA, scATAC, N, column_order = NULL, row_order = NULL) {
  
  # Get the unique clusters from both datasets
  scRNA_clusters <- unique(scRNA$cluster)
  scATAC_clusters <- unique(scATAC$cluster)
  
  # Initialize a matrix to store the p-values
  p_value_matrix <- matrix(NA, nrow = length(scRNA_clusters), ncol = length(scATAC_clusters),
                           dimnames = list(scRNA_clusters, scATAC_clusters))
  
  # Loop through each combination of scRNA and scATAC clusters
  for (scRNA_cluster in scRNA_clusters) {
    for (scATAC_cluster in scATAC_clusters) {
      
      # Get the marker genes for the current scRNA cluster
      scRNA_genes <- scRNA %>% filter(cluster == scRNA_cluster) %>% pull(gene)
      
      # Get the marker genes for the current scATAC cluster
      scATAC_genes <- scATAC %>% filter(cluster == scATAC_cluster) %>% pull(gene)
      
      # Calculate the values for the hypergeometric test
      K <- length(scRNA_genes)  # Number of marker genes in scRNA cluster
      n <- length(scATAC_genes) # Number of marker genes in scATAC cluster
      x <- length(intersect(scRNA_genes, scATAC_genes)) # Number of overlapping genes
      
      # Perform the hypergeometric test
      if (x > 0) {
        p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)
      } else {
        p_value <- 1 # If no overlap, p-value is 1
      }
      
      # Store the result in the matrix
      p_value_matrix[as.character(scRNA_cluster), as.character(scATAC_cluster)] <- p_value
    }
  }
  
  # Step 2: Transform the p-values to -log10(p-value) for better visualization
  log_p_value_matrix <- -log10(p_value_matrix)
  
  # Step 3: Reorder columns based on the specified order
  if (!is.null(column_order)) {
    log_p_value_matrix <- log_p_value_matrix[, column_order, drop = FALSE]
  }
  
  # Step 4: Reorder rows based on the specified order
  if (!is.null(row_order)) {
    log_p_value_matrix <- log_p_value_matrix[row_order, , drop = FALSE]
  }
  
  # Step 5: Create a heatmap with ComplexHeatmap using the log-transformed p-values
  heatmap <- Heatmap(log_p_value_matrix, name = "-log10(p-value)",
                     col = colorRamp2(c(0, max(log_p_value_matrix, na.rm = TRUE)), c("white", "red")),
                     cluster_rows = FALSE,  # Disable clustering of rows
                     cluster_columns = FALSE,  # Disable clustering of columns
                     row_order = row_order,  # Explicit row order
                     column_order = column_order,  # Explicit column order
                     row_names_gp = gpar(fontsize = 12), 
                     column_names_gp = gpar(fontsize = 12))
  
  # Step 6: Return a list with the matrix and heatmap
  return(list(
    p_value_matrix = p_value_matrix,
    heatmap = heatmap
  ))
}
