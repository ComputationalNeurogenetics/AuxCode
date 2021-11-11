# Some additional functions ----

# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  obj[['RNA_name']] <- obj[['RNA']]
  RNA <- obj@assays$RNA_name
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
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
    return(gene.part)
  } else {
    genes.part <- strsplit(gene.part,split = "::")[[1]]
     # Strip (ver something) away
    genes.part <- gsub(pattern = "\\(.+\\)", x = genes.part, replacement = "")
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
