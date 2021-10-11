# Some additional functions ----

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
