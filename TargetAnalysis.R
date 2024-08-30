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
    dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
    dbExecute(con, paste("CREATE TEMP TABLE filtered_tobias AS
    SELECT * FROM Tobias WHERE TF_gene_name = '",toupper(TF_name),"';",sep=""))
    
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
    li.end - li.start AS feature_center
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