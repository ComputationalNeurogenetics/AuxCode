extract_gene_enrichments <- function(ranks,gene_sets){
  # Plot enrichment and extract data
  enrichment_plot <- plotEnrichment(gene_sets, ranks) + 
    ggtitle("Enrichment Plot")
  
  # Extract the enrichment score profile data from the plot
  enrichment_data <- ggplot_build(enrichment_plot)$data[[1]]
  
  # Identify the leading and trailing edges
  positive_peak_index <- enrichment_data$x[which.max(enrichment_data$y)]
  negative_trough_index <- enrichment_data$x[which.min(enrichment_data$y)]
  
  # Ensure the trailing edge is calculated correctly
  trailing_edge_genes <- names(ranks)[(negative_trough_index + 1):length(ranks)]
  trailing_edge_genes <- trailing_edge_genes[trailing_edge_genes %in% gene_sets]
  
  # Leading edge genes
  leading_edge_genes <- names(ranks)[1:positive_peak_index]
  leading_edge_genes <- leading_edge_genes[leading_edge_genes %in% gene_sets]
  
  # Combine leading and trailing-edge genes
  combined_edges <- list(
    positive_leading_edge = leading_edge_genes,
    negative_trailing_edge = trailing_edge_genes
  )
  
  # Print combined leading and trailing-edge genes
  return(list(combined_edges=combined_edges, positive_peak_index=positive_peak_index, negative_peak_index=negative_trough_index))
}


drawTargetLinePlot <- function(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups){
  
  deg.data <- DEG.res$d1 %>% filter(p_val_adj<0.01) %>% mutate(rank = row_number()) %>% dplyr::select(ensg_id,rank,matches("exp_avg_log2$")) %>% pivot_longer(cols = where(is.double), names_to = "cell_group", values_to = "exp_avg_log2")
  
  limits.to.plot <- DEG.res$limits
  upper.limit=limits.to.plot$upper.limit
  lower.limit=limits.to.plot$lower.limit
  
  upper.pos=which.min(abs(filter(DEG.res$d1,p_val_adj<0.01) %>% pull(avg_log2FC) - upper.limit))
  lower.pos=which.min(abs(filter(DEG.res$d1,p_val_adj<0.01) %>% pull(avg_log2FC) - lower.limit))
  xmax.right=nrow(filter(DEG.res$d1,p_val_adj<0.01))
  
  cell_group_colors <- c("GA1_2_exp_avg_log2" = "blue", "GL1_2_exp_avg_log2" = "red") 
  
  l.1 <- ggplot(deg.data, aes(x=rank,y=exp_avg_log2, group=cell_group)) + geom_line(aes(colour = cell_group),alpha=.25) + theme_minimal() + theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 6),
    legend.position = "bottom") +
    scale_color_manual(values = cell_group_colors) +
    xlab("DEG target genes in order of avg_log2FC") + ylab("Expression avg log2") +
    geom_rect(aes(xmin=0,xmax=upper.pos,ymin=-0.1,ymax=0.1), fill="blue", alpha=0.002) +
    geom_rect(aes(xmin=lower.pos,xmax=xmax.right,ymin=-0.1,ymax=0.1), fill="red", alpha=0.002)

return(l.1)
}

drawTargetHeatmap <- function(DEG.res, target.genes, assay.averages,gene_name2id, gene_id2name,limits.to.plot, cell.groups){
  
  mid.upper.limit=limits.to.plot$mid.upper.limit
  mid.lower.limit=limits.to.plot$mid.lower.limit
  upper.limit=limits.to.plot$upper.limit
  lower.limit=limits.to.plot$lower.limit
  
  target.genes.in.order <- names(sort(DEG.res[target.genes], decreasing = TRUE))
  Targets.expression.avg <- log1p(t(assay.averages$RNA[sapply(target.genes.in.order,function(id){gene_name2id[[id]]}),cell.groups]))
  colnames(Targets.expression.avg) <- sapply(colnames(Targets.expression.avg),function(id){gene_id2name[[id]]})
  col_fun = colorRamp2(c(min(Targets.expression.avg), median(Targets.expression.avg),max(Targets.expression.avg)), c("white", "red","darkred"))
  
  color_factor <- factor(
    ifelse(DEG.res[target.genes.in.order] > upper.limit, "GA leading edge",
           ifelse(DEG.res[target.genes.in.order] <= mid.lower.limit & DEG.res[target.genes.in.order] >= mid.upper.limit, "Neutral",
                  ifelse(DEG.res[target.genes.in.order] < lower.limit, "GL trailing edge", "white")
           )
    ),
    levels = c("GL trailing edge", "white", "Neutral", "GA leading edge")
  )
  
  col_anno <- HeatmapAnnotation(target_gene_group = color_factor,
                                col = list(target_gene_group = c("GL trailing edge"="red" , "white" = "white", "Neutral" = "green", "GA leading edge" = "blue")),
                                annotation_legend_param = list(target_gene_group = list(title = "Target gene group")))
  
  h1 <- Heatmap(as.matrix(Targets.expression.avg), cluster_columns = FALSE,col=col_fun, show_column_names = FALSE,bottom_annotation = col_anno,heatmap_legend_param = list(title = "Target Expression (log1p)"))
  return(h1)
}

mean.fxn <- function (x,pseudocount.use=1, base=2)
{
  return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use,
             base = base))
}


getTargets_SR <- function(TF_name,con){
  if (TF_name!="Tal1"){
  dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
  dbExecute(con, paste("CREATE TEMP TABLE filtered_tobias AS
    SELECT * FROM Tobias WHERE TF_gene_name = '",toupper(TF_name),"';",sep=""))
  
  targetQuery <- paste("SELECT 
    gm.gene_name,
    gm.ensg_id,
    COUNT(CASE WHEN li.zscore > 0 THEN 1 ELSE NULL END) AS count_zscore_positive,
    COUNT(CASE WHEN li.zscore < 0 THEN 1 ELSE NULL END) AS count_zscore_negative,
    COUNT(DISTINCT li.feature) AS count_distinct_feature,
    COUNT(li.feature) AS count_feature,
    SUM(SR1_bound) as sum_SR1_bound,
    SUM(SR2_bound) as sum_SR2_bound,
    SUM(SR3_bound) as sum_SR3_bound,
    AVG(tb.SR1_score) AS SR1_score_average,
    AVG(tb.SR2_score) AS SR2_score_average,
    AVG(tb.SR3_score) AS SR3_score_average,
    AVG(tb.SR1_score) -  AVG(tb.SR2_score) AS SR1_SR2_increase,
    AVG(tb.SR2_score) -  AVG(tb.SR3_score) AS SR2_SR3_increase,
    AVG(tb.SR1_score) / NULLIF(AVG(tb.SR2_score), 0) AS SR1_SR_2_FT_ratio,
    AVG(tb.SR2_score) / NULLIF(AVG(tb.SR3_score), 0) AS SR2_SR_3_FT_ratio
FROM 
    links_s AS li
    JOIN filtered_tobias AS tb ON tb.features = li.feature
      JOIN (
        SELECT DISTINCT feature, target_gene_name
        FROM CT_data
        WHERE target_gene_name = '", TF_name, "'
      ) AS ct ON ct.feature = tb.features
    JOIN gene_metadata AS gm ON gm.ensg_id = li.ensg_id
WHERE
    tb.TF_gene_name = '",toupper(TF_name),"'
    AND ABS(li.zscore) > 2 
    AND li.pvalue < 0.01 
    AND (tb.SR1_bound = 1 OR tb.SR2_bound = 1 OR tb.SR3_bound = 1)
    AND tb.w_mean_cons > 0.5
GROUP BY 
    gm.gene_name;", sep="")
  } else if (TF_name=="Tal1") {
    dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
    dbExecute(con, "CREATE TEMP TABLE filtered_tobias AS
SELECT 
    *,
    (CASE 
        WHEN TFBS_name = 'TAL1.H12CORE.0.P.B_TAL1.H12CORE.0.P.B' THEN 1 
        ELSE 0 
     END) AS cond_0,
    (CASE 
        WHEN TFBS_name = 'TAL1.H12CORE.1.P.B_TAL1.H12CORE.1.P.B' THEN 1 
        ELSE 0 
     END) AS cond_1,
    (CASE 
        WHEN TFBS_name = 'TAL1.H12CORE.2.P.B_TAL1.H12CORE.2.P.B' THEN 1 
        ELSE 0 
     END) AS cond_2
FROM Tobias
WHERE TF_gene_name = 'TAL1'
  AND (TFBS_name IN ('TAL1.H12CORE.0.P.B_TAL1.H12CORE.0.P.B', 
                     'TAL1.H12CORE.1.P.B_TAL1.H12CORE.1.P.B', 
                     'TAL1.H12CORE.2.P.B_TAL1.H12CORE.2.P.B'));")
    
    targetQuery <- "SELECT 
    gm.gene_name,
    gm.ensg_id,
    COUNT(CASE WHEN li.zscore > 0 THEN 1 ELSE NULL END) AS count_zscore_positive,
    COUNT(CASE WHEN li.zscore < 0 THEN 1 ELSE NULL END) AS count_zscore_negative,
    COUNT(DISTINCT li.feature) AS count_distinct_feature,
    COUNT(li.feature) AS count_feature,
    SUM(SR1_bound) as sum_SR1_bound,
    SUM(SR2_bound) as sum_SR2_bound,
    SUM(SR3_bound) as sum_SR3_bound,
    AVG(tb.SR1_score) AS SR1_score_average,
    AVG(tb.SR2_score) AS SR2_score_average,
    AVG(tb.SR3_score) AS SR3_score_average,
    AVG(tb.SR1_score) -  AVG(tb.SR2_score) AS SR1_SR2_increase,
    AVG(tb.SR2_score) -  AVG(tb.SR3_score) AS SR2_SR3_increase,
    AVG(tb.SR1_score) / NULLIF(AVG(tb.SR2_score), 0) AS SR1_SR_2_FT_ratio,
    AVG(tb.SR2_score) / NULLIF(AVG(tb.SR3_score), 0) AS SR2_SR_3_FT_ratio,
    SUM(cond_0) AS TAL1_H12CORE_0_P_B_count,
    SUM(cond_1) AS TAL1_H12CORE_1_P_B_count,
    SUM(cond_2) AS TAL1_H12CORE_2_P_B_count
FROM 
    links_s AS li
    JOIN filtered_tobias AS tb ON tb.features = li.feature
      JOIN (
        SELECT DISTINCT feature, target_gene_name
        FROM CT_data
        WHERE target_gene_name = 'Tal1'
      ) AS ct ON ct.feature = tb.features
    JOIN gene_metadata AS gm ON gm.ensg_id = li.ensg_id
WHERE
    tb.TF_gene_name = 'TAL1'
    AND ABS(li.zscore) > 2 
    AND li.pvalue < 0.01 
    AND (tb.SR1_bound = 1 OR tb.SR2_bound = 1 OR tb.SR3_bound = 1)
    AND tb.w_mean_cons > 0.5
GROUP BY 
    gm.gene_name;"
  }
  
  target.data <- as_tibble(dbGetQuery(con, targetQuery))
  targets <- unique(pull(target.data, gene_name))
  targets_gene_id <- unique(pull(target.data, ensg_id))
  return(list(targets_gene_name=targets, targets_gene.id=targets_gene_id, target.data=target.data))
}

getTargets_rV2 <- function(TF_name,con){
  if (TF_name!="Tal1"){
    dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
    dbExecute(con, paste("CREATE TEMP TABLE filtered_tobias AS
    SELECT * FROM Tobias WHERE TF_gene_name = '",toupper(TF_name),"';",sep=""))
    
  targetQuery <- paste("SELECT 
    gm.gene_name,
    gm.ensg_id,
    COUNT(CASE WHEN li.zscore > 0 THEN 1 ELSE NULL END) AS count_zscore_positive,
    COUNT(CASE WHEN li.zscore < 0 THEN 1 ELSE NULL END) AS count_zscore_negative,
    COUNT(DISTINCT li.feature) AS count_distinct_feature,
    COUNT(li.feature) AS count_feature,
    SUM(GA1_2_bound) as sum_GA1_2_bound,
    SUM(GL1_2_bound) as sum_GL1_2_bound,
    AVG(tb.PRO1_2_score) AS PRO1_2_score_average,
    AVG(tb.CO1_2_score) AS CO1_2_score_average,
    AVG(tb.GA1_2_score) AS GA1_2_score_average,
    AVG(tb.GL1_2_score) AS GL1_2_score_average,
    AVG(tb.GA1_2_score) -  AVG(tb.PRO1_2_score) AS PRO1_2_GA_increase,
    AVG(tb.GL1_2_score) -  AVG(tb.PRO1_2_score) AS PRO1_2_GL_increase,
    AVG(tb.GA1_2_score) / NULLIF(AVG(tb.PRO1_2_score), 0) AS PRO1_2_GA_FT_ratio,
    AVG(tb.GL1_2_score) / NULLIF(AVG(tb.PRO1_2_score), 0) AS PRO1_2_GL_FT_ratio
FROM 
    links_s AS li
    JOIN filtered_tobias AS tb ON tb.features = li.feature
      JOIN (
        SELECT DISTINCT feature, target_gene_name
        FROM CT_data
        WHERE target_gene_name = '", TF_name, "'
      ) AS ct ON ct.feature = tb.features
    JOIN gene_metadata AS gm ON gm.ensg_id = li.ensg_id
WHERE
    tb.TF_gene_name = '",toupper(TF_name),"'
    AND ABS(li.zscore) > 2 
    AND li.pvalue < 0.01 
    AND (tb.GA1_2_bound = 1 OR tb.GL1_2_bound = 1)
    AND tb.w_mean_cons > 0.5
GROUP BY 
    gm.gene_name;", sep="")
  }
  else if (TF_name=="Tal1") {
    dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
    dbExecute(con, "CREATE TEMP TABLE filtered_tobias AS
SELECT 
    *,
    (CASE 
        WHEN TFBS_name = 'TAL1.H12CORE.0.P.B_TAL1.H12CORE.0.P.B' THEN 1 
        ELSE 0 
     END) AS cond_0,
    (CASE 
        WHEN TFBS_name = 'TAL1.H12CORE.1.P.B_TAL1.H12CORE.1.P.B' THEN 1 
        ELSE 0 
     END) AS cond_1,
    (CASE 
        WHEN TFBS_name = 'TAL1.H12CORE.2.P.B_TAL1.H12CORE.2.P.B' THEN 1 
        ELSE 0 
     END) AS cond_2
FROM Tobias
WHERE TF_gene_name = 'TAL1'
  AND (TFBS_name IN ('TAL1.H12CORE.0.P.B_TAL1.H12CORE.0.P.B', 
                     'TAL1.H12CORE.1.P.B_TAL1.H12CORE.1.P.B', 
                     'TAL1.H12CORE.2.P.B_TAL1.H12CORE.2.P.B'));")
    
    targetQuery <- "SELECT 
    gm.gene_name,
    gm.ensg_id,
    COUNT(CASE WHEN li.zscore > 0 THEN 1 ELSE NULL END) AS count_zscore_positive,
    COUNT(CASE WHEN li.zscore < 0 THEN 1 ELSE NULL END) AS count_zscore_negative,
    COUNT(DISTINCT li.feature) AS count_distinct_feature,
    COUNT(li.feature) AS count_feature,
    SUM(GA1_2_bound) as sum_GA1_2_bound,
    SUM(GL1_2_bound) as sum_GL1_2_bound,
    AVG(tb.PRO1_2_score) AS PRO1_2_score_average,
    AVG(tb.CO1_2_score) AS CO1_2_score_average,
    AVG(tb.GA1_2_score) AS GA1_2_score_average,
    AVG(tb.GL1_2_score) AS GL1_2_score_average,
    AVG(tb.GA1_2_score) -  AVG(tb.PRO1_2_score) AS PRO1_2_GA_increase,
    AVG(tb.GL1_2_score) -  AVG(tb.PRO1_2_score) AS PRO1_2_GL_increase,
    AVG(tb.GA1_2_score) / NULLIF(AVG(tb.PRO1_2_score), 0) AS PRO1_2_GA_FT_ratio,
    AVG(tb.GL1_2_score) / NULLIF(AVG(tb.PRO1_2_score), 0) AS PRO1_2_GL_FT_ratio,
    SUM(cond_0) AS TAL1_H12CORE_0_P_B_count,
    SUM(cond_1) AS TAL1_H12CORE_1_P_B_count,
    SUM(cond_2) AS TAL1_H12CORE_2_P_B_count
FROM 
    links_s AS li
    JOIN filtered_tobias AS tb ON tb.features = li.feature
    JOIN (
        SELECT DISTINCT feature, target_gene_name
        FROM CT_data
        WHERE target_gene_name = 'Tal1'
    ) AS ct ON ct.feature = tb.features
    JOIN gene_metadata AS gm ON gm.ensg_id = li.ensg_id
WHERE
    ABS(li.zscore) > 2 
    AND li.pvalue < 0.01 
    AND (tb.GA1_2_bound = 1 OR tb.GL1_2_bound = 1)
    AND tb.w_mean_cons > 0.5
GROUP BY 
    gm.gene_name;"
  }
  
  # AND ct.target_gene_name = '",TF_name,"'
  target.data <- as_tibble(dbGetQuery(con, targetQuery))
  targets <- unique(pull(target.data, gene_name))
  targets_gene_id <- unique(pull(target.data, ensg_id))
  return(list(targets_gene_name=targets, targets_gene.id=targets_gene_id, target.data=target.data))
}

getTargets_rV2_feature_level <- function(TF_name,con){
    TF_name.tobias <- TF_name
    if (TF_name=="Ebf1"){TF_name.tobias<-"COE1"}
    dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
    dbExecute(con, paste("CREATE TEMP TABLE filtered_tobias AS
    SELECT * FROM Tobias WHERE TF_gene_name = '",toupper(TF_name.tobias),"';",sep=""))
    
    targetQuery <- paste("SELECT 
    gm.gene_name,
    gm.ensg_id,
    tb.TFBS_name,
    li.feature AS feature,
    GA1_2_bound as GA1_2_bound,
    GL1_2_bound as GL1_2_bound,
    tb.PRO1_2_score,
    tb.CO1_2_score,
    tb.GA1_2_score,
    tb.GL1_2_score,
    tb.GA1_2_score -  tb.PRO1_2_score AS PRO1_2_GA_increase,
    tb.GL1_2_score -  tb.PRO1_2_score AS PRO1_2_GL_increase,
    tb.GA1_2_score / NULLIF(tb.PRO1_2_score, 0) AS PRO1_2_GA_FT_ratio,
    tb.GL1_2_score / NULLIF(tb.PRO1_2_score, 0) AS PRO1_2_GL_FT_ratio,
    li.zscore,
    ABS(((li.end+li.start)/2)-gm.TSS_start) AS linkage_dist
    
FROM 
    links_s AS li
    JOIN filtered_tobias AS tb ON tb.features = li.feature
      JOIN (
        SELECT DISTINCT feature, target_gene_name
        FROM CT_data
        WHERE target_gene_name = '", TF_name, "'
      ) AS ct ON ct.feature = tb.features
    JOIN gene_metadata AS gm ON gm.ensg_id = li.ensg_id
WHERE
    tb.TF_gene_name = '",toupper(TF_name.tobias),"'
    AND ABS(li.zscore) > 2 
    AND li.pvalue < 0.01 
    AND (tb.GA1_2_bound = 1 OR tb.GL1_2_bound = 1)
    AND tb.w_mean_cons > 0.5", sep="")

  target.data <- as_tibble(dbGetQuery(con, targetQuery))
  targets <- unique(pull(target.data, gene_name))
  targets_gene_id <- unique(pull(target.data, ensg_id))
  return(list(targets_gene_name=targets, targets_gene.id=targets_gene_id, target.data=target.data))
}

do.GSEA <- function(targets,DEG.res,TF.name, comp.title){
  # lower.thr <- .25
  # upper.thr <- .75
  
  DEG.res <- arrange(DEG.res, desc(avg_log2FC))
  DEG.res.GSEA <- DEG.res$avg_log2FC
  names(DEG.res.GSEA) <- DEG.res$gene_name
  
  edges <- extract_gene_enrichments(ranks=sort(DEG.res.GSEA,decreasing=TRUE), gene_sets = targets)
  
  # upper.limit <- quantile(DEG.res$avg_log2FC,upper.thr)
  # upper.pos <- which.min(abs(DEG.res$avg_log2FC-upper.limit))
  upper.pos <- edges$positive_peak_index
  lower.pos <- edges$negative_peak_index
  
  upper.limit <- DEG.res$avg_log2FC[upper.pos]
  lower.limit <- DEG.res$avg_log2FC[lower.pos]
  
  # lower.limit <- quantile(DEG.res$avg_log2FC,lower.thr)
  # lower.pos <- which.min(abs(DEG.res$avg_log2FC-lower.limit))
  # 
  mid.thr.width <- .25
  width.around.zero <- round(nrow(DEG.res)*(mid.thr.width/2), digits = 0)
  zero.index <- which.min(abs(DEG.res$avg_log2FC))
  # 
  mid.lower.index <- zero.index-width.around.zero
  mid.upper.index <- zero.index+width.around.zero
  mid.lower.limit <- DEG.res$avg_log2FC[mid.lower.index]
  mid.upper.limit <- DEG.res$avg_log2FC[mid.upper.index]
  #   
  xmax.right <- nrow(DEG.res)
  # 
  limits <- list(upper.pos=upper.pos,lower.pos=lower.pos,upper.limit=upper.limit,lower.limit=lower.limit,mid.lower.limit=mid.lower.limit,mid.upper.limit=mid.upper.limit)
  
  
  p1 <- plotEnrichment(targets, DEG.res.GSEA) + ggtitle(paste(TF.name,"target enrichments in", comp.title)) + geom_rect(aes(xmin=0,xmax=upper.pos,ymin=-0.05,ymax=0.05), fill="blue", alpha=0.002) + geom_rect(aes(xmin=lower.pos,xmax=xmax.right,ymin=-0.05,ymax=0.05), fill="red", alpha=0.002) #geom_rect(aes(xmin=mid.lower.index,xmax=mid.upper.index,ymin=-0.05,ymax=0.05), fill="green", alpha=0.002) + 
  d1 <- DEG.res %>% dplyr::filter(gene_name %in% targets)
  
  d1$target_gene_group <- factor(
    ifelse(d1$avg_log2FC > upper.limit, "blue",
           ifelse(d1$avg_log2FC >= mid.lower.limit & d1$avg_log2FC <= mid.upper.limit, "white",
                  ifelse(d1$avg_log2FC < lower.limit, "red", "white")
           )
    ),
    levels = c("red", "white", "blue")
  )
  
  return(list(p1=p1,d1=d1,limits=limits,edges=edges))
}
# Define the function
plot_heatmap_with_variability <- function(
    gene_list,            # List of gene names (e.g., Tal1.causal.targets)
    combined_tibble,      # Dataframe with the target information
    rna_data,             # Matrix of RNA expression data
    gene_id2name,         # A named vector to map gene IDs to gene names
    TF,                   # Required argument: TF to filter on (e.g., "Tal1")
    zscore_condition = "> 0", # User-defined zscore condition (e.g., "> 0" or "< 0")
    custom_filter = NULL, # Optional character string for additional filter (e.g., "& Gata2_ct_overlap == 0 & Gata3_ct_overlap == 0")
    filter_threshold = 1.2, # RNA expression filter threshold
    output_file = "heatmap_with_variability.pdf", # Output file for the heatmap
    column_samples = c("PRO1","PRO2","CO1","CO2","GA1","GA2","GA3","GA4","GA5","GA6","GL1","GL2","GL3","GL4","GL5"), # Column names (samples) to plot
    bin_ranges = c("0-5 kbp", "> 5 kbp - 50 kbp", "> 50 kbp") # Bin ranges for splitting rows
) {
  
  # Construct the filter expression dynamically, incorporating the zscore_condition
  filter_expression <- paste0("TF == '", TF, "' & zscore ", zscore_condition)
  
  # Add the custom filter string if provided
  if (!is.null(custom_filter)) {
    filter_expression <- paste0(filter_expression, " ", custom_filter)
  }
  
  # Filter the combined_tibble using the dynamically constructed filter expression
  pos.targets <- combined_tibble %>%
    filter(eval(parse(text = filter_expression))) # Dynamically evaluates the filtering
  
  # Group and mutate for bin_category
  pos.targets <- pos.targets %>%
    group_by(ensg_id) %>%
    summarize(linkage_kbp = min(linkage_kbp)) %>%
    mutate(bin_category = case_when(
      linkage_kbp <= 5 ~ bin_ranges[1],
      linkage_kbp > 5 & linkage_kbp <= 50 ~ bin_ranges[2],
      linkage_kbp > 50 ~ bin_ranges[3]
    )) %>%
    ungroup()
  
  # Filter expressed genes for heatmap
  exp.filt <- rownames(rna_data[rowSums(rna_data > filter_threshold) > 0,])
  pos.targets <- pos.targets %>% filter(ensg_id %in% exp.filt)
  
  # Data to plot (log transformed RNA expression)
  data2plot <- log1p(rna_data[pos.targets$ensg_id, column_samples])
  
  # Set gene names as row names
  rownames(data2plot) <- sapply(rownames(data2plot), function(g) { gene_id2name[[g]] })
  
  # Calculate variability (standard deviation) per gene
  gene_variability <- apply(data2plot, 1, sd)
  
  # Ensure the row order of data2plot and pos.targets matches, add variability to pos.targets
  pos.targets <- pos.targets %>%
    mutate(gene_name = sapply(ensg_id, function(g) { gene_id2name[[g]] })) %>%
    filter(gene_name %in% rownames(data2plot)) %>%
    arrange(factor(gene_name, levels = rownames(data2plot))) %>%
    mutate(variability = gene_variability[match(gene_name, rownames(data2plot))])
  
  # Ensure that gene_variability matches the row names of data2plot
  gene_variability <- gene_variability[match(rownames(data2plot), names(gene_variability))]
  
  # Calculate average variability per bin_category
  avg_variability_by_bin <- pos.targets %>%
    group_by(bin_category) %>%
    summarise(avg_variability = mean(variability, na.rm = TRUE))
  
  # Append average variability to bin_category names
  pos.targets <- pos.targets %>%
    left_join(avg_variability_by_bin, by = "bin_category") %>%
    mutate(bin_category_label = paste0(bin_category, " (Avg Var: ", round(avg_variability, 3), ")"))
  
  # Create row annotation (binary and variability)
  row_annotation <- rowAnnotation(
    Causal_Target = anno_simple(as.integer(rownames(data2plot) %in% gene_list), col = c("0" = "white", "1" = "green")),
    Variability = anno_simple(gene_variability, col = colorRamp2(c(min(gene_variability), max(gene_variability)), c("blue", "red")))
  )
  
  # Color scale for the heatmap
  color_scale <- colorRamp2(c(min(data2plot), max(data2plot)), c("white", "red"))
  
  # Plot the heatmap
  pdf(file = output_file, width = 10, height = 25)
  h1 <- Heatmap(
    data2plot,
    name = "RNA expression",
    row_split = pos.targets$bin_category_label,  # Use modified bin_category labels with variability
    cluster_rows = TRUE,                         # Cluster rows within each bin_category
    cluster_columns = FALSE,                     # Remove column clustering
    show_row_names = TRUE,                       # Show row names
    show_column_names = TRUE,                    # Show column names
    clustering_distance_rows = "euclidean",
    row_names_gp = gpar(fontsize = 6),
    clustering_method_rows = "ward.D2",
    col = color_scale,
    left_annotation = row_annotation              # Add row annotation to the heatmap
  )
  print(h1)
  dev.off()
}


# Step 1: Adding Transfer Entropy Results to the input_tibble
add_TE_to_tibble <- function(input_tibble, te_results_list) {
  # Create a dataframe to hold TE results
  te_df <- do.call(rbind, lapply(names(te_results_list), function(gene) {
    te_res <- te_results_list[[gene]]$coef
    data.frame(
      ensg_id = gene,
      X_Y_te_p_value = te_res["X->Y", "p-value"],
      Y_X_te_p_value = te_res["Y->X", "p-value"],
      X_Y_te_coef = te_res["X->Y", "te"],
      Y_X_te_coef = te_res["Y->X", "te"],
      stringsAsFactors = FALSE
    )
  }))
  
  # Merge TE results with the original tibble by ensg_id
  input_tibble_with_te <- input_tibble %>%
    left_join(te_df, by = "ensg_id")
  
  return(input_tibble_with_te)
}

# Step 1: Adding Transfer Entropy Results to the input_tibble
add_TE_to_tibble <- function(input_tibble, te_results_list) {
  # Create a dataframe to hold TE results
  te_df <- do.call(rbind, lapply(names(te_results_list), function(gene) {
    te_res <- te_results_list[[gene]]$coef
    data.frame(
      ensg_id = gene,
      X_Y_te_p_value = te_res["X->Y", "p-value"],
      Y_X_te_p_value = te_res["Y->X", "p-value"],
      X_Y_te_coef = te_res["X->Y", "te"],
      Y_X_te_coef = te_res["Y->X", "te"],
      stringsAsFactors = FALSE
    )
  }))
  
  # Merge TE results with the original tibble by ensg_id
  input_tibble_with_te <- input_tibble %>%
    left_join(te_df, by = "ensg_id")
  
  return(input_tibble_with_te)
}

# Define the function to process TF targets with optimized thresholds
process_TF_targets_with_TE <- function(input_tibble_with_te, TF) {
  
  # Filter based on the specified TF
  filtered_targets <- input_tibble_with_te %>%
    filter(TF == !!TF)
  
  # Classify genes based on TE results, correlations, and log2FC
  filtered_targets <- filtered_targets %>%
    mutate(
      regulation = case_when(
        zscore > 0 & avg_log2FC > 0.1 ~ "Activation",  # TF likely activates this gene
        zscore < 0 & avg_log2FC < -0.1 ~ "Repression", # TF likely represses this gene
        TRUE ~ "Unknown"                               # Mixed or weak signal
      ),
      causality_present_X_Y = ifelse(X_Y_te_p_value < 0.05 & X_Y_te_coef > 0.02, TRUE, FALSE), # X->Y TE significant and strong
      causality_present_Y_X = ifelse(Y_X_te_p_value < 0.05 & Y_X_te_coef > 0.02, TRUE, FALSE), # Y->X TE significant and strong
      TE_relationship = case_when(
        causality_present_X_Y & !causality_present_Y_X ~ "TF->Target",
        causality_present_Y_X & !causality_present_X_Y ~ "Target->TF",
        causality_present_X_Y & causality_present_Y_X ~ "Bidirectional",
        TRUE ~ "No significant TE"
      ),
      # Correlation thresholds for PRO-GA and PRO-GL transitions
      positive_cor_GA = ifelse(TF.cor.PRO.GA > 0.3, TRUE, FALSE),
      negative_cor_GA = ifelse(TF.cor.PRO.GA < -0.3, TRUE, FALSE),
      positive_cor_GL = ifelse(TF.cor.PRO.GL > 0.3, TRUE, FALSE),
      negative_cor_GL = ifelse(TF.cor.PRO.GL < -0.3, TRUE, FALSE),
      has_Gata2_overlap = ifelse(TF == "Tal1" & Gata2_ct_overlap > 0, TRUE, FALSE), # Tal1-specific Gata2 overlap
      has_Gata3_overlap = ifelse(TF == "Tal1" & Gata3_ct_overlap > 0, TRUE, FALSE), # Tal1-specific Gata3 overlap
      proximal_binding = ifelse(linkage_kbp <= 50, TRUE, FALSE),                    # Binding proximity < 50kb
      expressed_in_GA = ifelse(pct.1 > 0.1, TRUE, FALSE),                           # Expression in GA > 10%
      expressed_in_GL = ifelse(pct.2 > 0.1, TRUE, FALSE),                           # Expression in GL > 10%
      PRO_GA_increase = ifelse(PRO1_2_GA_FT_ratio > 1, TRUE, FALSE),                # Footprint increase PRO->GA
      PRO_GL_increase = ifelse(PRO1_2_GL_FT_ratio > 1, TRUE, FALSE)                 # Footprint increase PRO->GL
    )
  
  # Query function that incorporates thresholds for TE, correlation, and log2FC
  query_TF_targets <- function(type = c("Activation", "Repression"),
                               te_direction = NULL, condition = NULL,
                               Gata2_overlap = NULL, Gata3_overlap = NULL,
                               proximal_binding = NULL, PRO_GA_increase = NULL,
                               PRO_GL_increase = NULL, expressed_in_GA = NULL,
                               expressed_in_GL = NULL, strong_TE = NULL, lineage_focus = NULL,
                               positive_cor_GA = NULL, negative_cor_GA = NULL,
                               positive_cor_GL = NULL, negative_cor_GL = NULL) {
    
    # Start with basic filtering on activation/repression
    query_result <- filtered_targets %>%
      filter(regulation %in% type)
    
    # Add filtering by TE direction if specified
    if (!is.null(te_direction)) {
      query_result <- query_result %>%
        filter(TE_relationship == te_direction)
    }
    
    # Filter by condition-specific regulation if specified
    if (!is.null(condition)) {
      query_result <- query_result %>%
        filter(condition_specificity == condition)
    }
    
    # Add filtering for Gata2 overlap if specified and the TF is Tal1
    if (!is.null(Gata2_overlap) & TF == "Tal1") {
      query_result <- query_result %>%
        filter(has_Gata2_overlap == Gata2_overlap)
    }
    
    # Add filtering for Gata3 overlap if specified and the TF is Tal1
    if (!is.null(Gata3_overlap) & TF == "Tal1") {
      query_result <- query_result %>%
        filter(has_Gata3_overlap == Gata3_overlap)
    }
    
    # Filter by proximity to binding site if specified
    if (!is.null(proximal_binding)) {
      query_result <- query_result %>%
        filter(proximal_binding == proximal_binding)
    }
    
    # Filter by footprint score increase from PRO to GA or GL
    if (!is.null(PRO_GA_increase)) {
      query_result <- query_result %>%
        filter(PRO_GA_increase == PRO_GA_increase)
    }
    if (!is.null(PRO_GL_increase)) {
      query_result <- query_result %>%
        filter(PRO_GL_increase == PRO_GL_increase)
    }
    
    # Filter by gene expression in GA or GL if specified
    if (!is.null(expressed_in_GA)) {
      query_result <- query_result %>%
        filter(expressed_in_GA == expressed_in_GA)
    }
    if (!is.null(expressed_in_GL)) {
      query_result <- query_result %>%
        filter(expressed_in_GL == expressed_in_GL)
    }
    
    # Filter by strong TE coefficient if specified
    if (!is.null(strong_TE)) {
      if (strong_TE == "X_Y") {
        query_result <- query_result %>%
          filter(causality_present_X_Y == TRUE)
      } else if (strong_TE == "Y_X") {
        query_result <- query_result %>%
          filter(causality_present_Y_X == TRUE)
      }
    }
    
    # Filter by lineage focus (GA or GL)
    if (!is.null(lineage_focus)) {
      if (lineage_focus == "GA") {
        query_result <- query_result %>%
          filter(PRO_GA_increase == TRUE)
      } else if (lineage_focus == "GL") {
        query_result <- query_result %>%
          filter(PRO_GL_increase == TRUE)
      }
    }
    
    # Add filtering for positive or negative correlations in GA or GL
    if (!is.null(positive_cor_GA)) {
      query_result <- query_result %>%
        filter(positive_cor_GA == positive_cor_GA)
    }
    if (!is.null(negative_cor_GA)) {
      query_result <- query_result %>%
        filter(negative_cor_GA == negative_cor_GA)
    }
    if (!is.null(positive_cor_GL)) {
      query_result <- query_result %>%
        filter(positive_cor_GL == positive_cor_GL)
    }
    if (!is.null(negative_cor_GL)) {
      query_result <- query_result %>%
        filter(negative_cor_GL == negative_cor_GL)
    }
    
    return(query_result)
  }
  
  # Return both the processed data and the query function
  list(
    data = filtered_targets,
    query_function = query_TF_targets
  )
}

