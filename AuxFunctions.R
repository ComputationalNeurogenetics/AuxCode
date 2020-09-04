# Some additional functions ----

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
  TF.names.from.enrichments <- pull(enriched.tb.filt,Term)
  
  # Intersect names ----
  intersected.TF.names <- intersect(TF.names.from.enrichments,TF.names.from.motifs)
  
  # Pull out genes present both in markers and enriched among the downstream target genes of each intersected TF ----
  program.genes <- lapply(intersected.TF.names, function(term){
    # Pulling out genes of each term from the intersect
    genes.downstream <- strsplit(x = pull(filter(enriched.tb.filt, Term==term), Genes),split=";")[[1]] 
    # Intersecting genes downstream with markers
    genes.for.program <- intersect(genes.downstream,toupper(markers))
    # Add TF itself to the program
    str_to_title(c(genes.for.program, term))
  })
  # Format output ----
  names(program.genes) <- intersected.TF.names
  return(list(programs=program.genes, db.name=db.name, adj.p.value.threshold = adj.p.value))
}

Fetch_Enrichments <- function(genes,db.name){
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE else websiteLive <- TRUE
  if (websiteLive) {
    enriched <- enrichr(genes, db.name)
    enriched <- lapply(enriched, as_tibble)
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