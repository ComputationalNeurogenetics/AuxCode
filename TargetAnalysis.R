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
  
  # Prepare the data
  deg.data <- DEG.res$d1 %>%
    filter(p_val_adj < 0.01) %>%
    mutate(rank = row_number()) %>%
    dplyr::select(ensg_id, rank, matches("exp_avg_log2$")) %>%
    pivot_longer(cols = where(is.double), names_to = "cell_group", values_to = "exp_avg_log2")
  
  # Extract limits for plotting
  upper.limit = limits.to.plot$upper.limit
  lower.limit = limits.to.plot$lower.limit
  
  # Find positions for upper and lower limits
  upper.pos = which.min(abs(filter(DEG.res$d1, p_val_adj < 0.01) %>% pull(avg_log2FC) - upper.limit))
  lower.pos = which.min(abs(filter(DEG.res$d1, p_val_adj < 0.01) %>% pull(avg_log2FC) - lower.limit))
  xmax.right = nrow(filter(DEG.res$d1, p_val_adj < 0.01))
  
  # Define colors for cell groups
  cell_group_colors <- c("GA1_2_exp_avg_log2" = "blue", "GL1_2_exp_avg_log2" = "red")
  
  # Create the line plot with limits applied
  l.1 <- ggplot(deg.data, aes(x = rank, y = exp_avg_log2, group = cell_group)) +
    geom_line(aes(colour = cell_group), alpha = 0.55) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 6),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = cell_group_colors) +
    xlab("DEG target genes in order of avg_log2FC") +
    ylab("Expression avg log2") +
    # Fixing the geom_rect aesthetics by using fixed values
    geom_rect(xmin = 0, xmax = upper.pos, ymin = -0.1, ymax = 0.1, fill = "blue", alpha = 0.002) +
    geom_rect(xmin = lower.pos, xmax = xmax.right, ymin = -0.1, ymax = 0.1, fill = "red", alpha = 0.002)
  
  return(l.1)
}

library(ggplot2)
library(dplyr)

drawEdgeLinePlot <- function(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups, edge_type = c("leading", "trailing")){
  
  # Match the edge_type argument to ensure valid input
  edge_type <- match.arg(edge_type)
  
  # Filter for leading or trailing edge genes
  if (edge_type == "leading") {
    target_genes <- DEG.res$edges$combined_edges$positive_leading_edge
    color_fill <- "blue"
  } else if (edge_type == "trailing") {
    target_genes <- DEG.res$edges$combined_edges$negative_trailing_edge
    color_fill <- "red"
  }
  
  # Filter the data for selected edge type genes
  deg.data <- DEG.res$d1 %>%
    filter(gene_name %in% target_genes) %>%
    mutate(rank = row_number()) %>%
    dplyr::select(ensg_id, rank, matches("exp_avg_log2$")) %>%
    pivot_longer(cols = where(is.double), names_to = "cell_group", values_to = "exp_avg_log2")
  
  # Extract limits for plotting
  upper.limit = limits.to.plot$upper.limit
  lower.limit = limits.to.plot$lower.limit
  
  xmax.right = nrow(deg.data)
  
  # Define colors for cell groups
  cell_group_colors <- c("GA1_2_exp_avg_log2" = "blue", "GL1_2_exp_avg_log2" = "red")
  
  # Create the line plot with limits applied
  l.1 <- ggplot(deg.data, aes(x = rank, y = exp_avg_log2, group = cell_group)) +
    geom_line(aes(colour = cell_group), alpha = 0.55) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 6),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = cell_group_colors) +
    xlab(paste(edge_type, "edge target genes in order of avg_log2FC")) +
    ylab("Expression avg log2") +
    # Add rectangular highlighting if needed, can be removed if not required
    geom_rect(xmin = 0, xmax = xmax.right, ymin = -0.1, ymax = 0.1, fill = color_fill, alpha = 0.002)
  
  return(l.1)
}

plotEdgeScatterWithRawValues <- function(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups, edge_type = c("leading", "trailing"), log_scale = FALSE) {
  
  # Match the edge_type argument to ensure valid input
  edge_type <- match.arg(edge_type)
  
  # Filter for leading or trailing edge genes
  if (edge_type == "leading") {
    target_genes <- DEG.res$edges$combined_edges$positive_leading_edge
    color_fill <- "blue"
  } else if (edge_type == "trailing") {
    target_genes <- DEG.res$edges$combined_edges$negative_trailing_edge
    color_fill <- "red"
  }
  
  # Filter the DEG data for the selected edge type genes
  deg.data <- DEG.res$d1 %>%
    filter(gene_name %in% target_genes) %>%
    mutate(rank = row_number()) %>%
    dplyr::select(ensg_id, gene_name, rank, avg_log2FC, GA1_2_exp_avg, GL1_2_exp_avg)  # Use raw expression values
  
  # Apply log scale if needed
  if (log_scale) {
    deg.data <- deg.data %>%
      mutate(GA1_2_exp_avg = log2(GA1_2_exp_avg + 1),
             GL1_2_exp_avg = log2(GL1_2_exp_avg + 1))
  }
  
  # Extract the limits for plotting
  upper.limit <- limits.to.plot$upper.limit
  lower.limit <- limits.to.plot$lower.limit
  mid.upper.limit <- limits.to.plot$mid.upper.limit
  mid.lower.limit <- limits.to.plot$mid.lower.limit
  
  # Create the scatter plot with limits applied
  l1 <- ggplot(deg.data, aes(x = GL1_2_exp_avg, y = GA1_2_exp_avg)) +
    geom_point(color = color_fill, size = 2, alpha = 0.6) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.position = "none"
    ) +
    xlab("GL1_2 Expression log2(avg+1)") +
    ylab("GA1_2 Expression log2(avg+1)") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = c(upper.limit, mid.upper.limit, mid.lower.limit, lower.limit), linetype = "dotted", color = "gray") +
    geom_vline(xintercept = c(upper.limit, mid.upper.limit, mid.lower.limit, lower.limit), linetype = "dotted", color = "gray") +
    annotate("rect", xmin = lower.limit, xmax = upper.limit, ymin = lower.limit, ymax = upper.limit, fill = color_fill, alpha = 0.05) +
    geom_text_repel(
      aes(label = gene_name),
      max.overlaps = 30,  # Limits number of labels
      force = 10,         # Adjusts how strong the repelling force is
      box.padding = 0.15,  # Adjusts spacing around labels
      point.padding = 0.5,
      segment.color = "grey50",  # Color of the connecting lines
      segment.size = 0.25,        # Thickness of the connecting lines
      segment.alpha=0.75,
      size = 3,
      max.time = 1,       # Time for placing labels
      seed = 42           # Ensures reproducibility
    )
  
  return(l1)
}

drawTargetBarPlot <- function(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups, limit_x_axis = TRUE) {
  
  # Filter the significant genes based on p_val_adj < 0.01
  deg.data <- DEG.res$d1 %>%
    filter(p_val_adj < 0.01) %>%
    dplyr::select(ensg_id, avg_log2FC, matches("exp_avg_log2$"))
  
  # Apply the x-axis limits based on avg_log2FC
  if (limit_x_axis) {
    deg.data <- deg.data %>%
      filter(avg_log2FC >= limits.to.plot$lower.limit & avg_log2FC <= limits.to.plot$upper.limit)
  }
  
  # Add rank based on the filtered data (after applying limits)
  deg.data <- deg.data %>%
    mutate(rank = row_number()) %>%
    pivot_longer(cols = where(is.double), names_to = "cell_group", values_to = "exp_avg_log2") 
  
  # Define colors for cell groups
  cell_group_colors <- c("GA1_2_exp_avg_log2" = "blue", "GL1_2_exp_avg_log2" = "red")
  
  # Create the stacked bar plot
  l.1 <- ggplot(deg.data, aes(x = rank, y = exp_avg_log2, fill = cell_group)) +
    geom_bar(stat = "identity", position = "stack", alpha = .75) +  # Stacked bar plot
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 6),
      legend.position = "bottom"
    ) +
    scale_fill_manual(values = cell_group_colors) +
    xlab("DEG target genes in order of avg_log2FC") +
    ylab("Expression avg log2")
  
  return(l.1)
}

# Function to create a volcano plot with limits based on drawTargetBarPlot limits.to.plot
drawVolcanoPlotWithLimits <- function(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups) {
  
  # Filter significant genes based on p_val_adj < 0.01
  deg.data <- DEG.res$d1 %>%
    filter(p_val_adj < 0.01) %>%
    arrange(desc(avg_log2FC))  # Arrange by avg_log2FC
  
  # Extract limits for plotting
  upper.limit <- limits.to.plot$upper.limit
  lower.limit <- limits.to.plot$lower.limit
  
  # Find the index positions for the limits in the sorted data
  upper.pos <- which.min(abs(deg.data$avg_log2FC - upper.limit))  # Last gene in the blue region
  lower.pos <- which.min(abs(deg.data$avg_log2FC - lower.limit))  # First gene in the red region
  
  # Select the genes that are outside of the defined regions
  selected_genes <- deg.data %>%
    filter(row_number() <= upper.pos | row_number() >= lower.pos)
  
  # Create the volcano plot
  volcano_plot <- ggplot(deg.data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(alpha = 0.5, color = "gray") +  # Background genes
    geom_point(data = selected_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC), size = 3) +  # Highlighted points
    scale_color_gradient(low = "red", high = "blue") +  # Color gradient for log fold change
    labs(
      x = "Log2 Fold Change",
      y = "-log10 Adjusted P-value",
      title = "Volcano Plot Highlighting Genes Outside Defined Limits"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  return(volcano_plot)
}

library(ggplot2)
library(dplyr)

# Function to draw scatterplot of GA expression vs GL expression for target genes, limited to blue and red regions
plotGA_vs_GL_ScatterWithLimits <- function(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups) {
  
  # Extract limits for plotting
  upper.limit <- limits.to.plot$upper.limit  # GA trailing edge limit (should be blue)
  lower.limit <- limits.to.plot$lower.limit  # GL trailing edge limit (should be red)
  
  # Filter the data for significant target genes (p_val_adj < 0.01)
  deg.data <- DEG.res$d1 %>%
    filter(p_val_adj < 0.01)
  
  # Separate GA trailing edge (blue) and GL trailing edge (red) region genes
  ga_trailing_edge_genes <- deg.data %>%
    filter(avg_log2FC <= upper.limit)  # GA trailing edge (blue)
  
  gl_trailing_edge_genes <- deg.data %>%
    filter(avg_log2FC >= lower.limit)  # GL trailing edge (red)
  
  # Combine GA and GL trailing edge genes and add a color label
  combined_genes <- bind_rows(
    ga_trailing_edge_genes %>% mutate(region = "GA trailing edge"),  # GA trailing edge (blue)
    gl_trailing_edge_genes %>% mutate(region = "GL trailing edge")   # GL trailing edge (red)
  )
  
  # Create the scatter plot comparing GL vs GA expression, colored by region
  scatter_plot <- ggplot(combined_genes, aes(x = GL1_2_exp_avg_log2, y = GA1_2_exp_avg_log2, color = region)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("GA trailing edge" = "blue", "GL trailing edge" = "red")) +  # Corrected region color mapping
    labs(
      x = "GL Expression (avg log2)",
      y = "GA Expression (avg log2)",
      title = "Scatter Plot: GA vs GL Expression (Trailing Edges)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(scatter_plot)
}

# Example usage:
# plotGA_vs_GL_ScatterWithLimits(DEG.res, gene_name2id, gene_id2name, limits.to.plot, cell.groups)


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

getTargets_rV2_feature_level_v2 <- function(TF_name, con, include_extra_columns = FALSE) {
  # Adjust TF name for Tobias DB if necessary
  TF_name.tobias <- TF_name
  if (TF_name == "Ebf1") {
    TF_name.tobias <- "COE1"
  }
  
  # Create temporary table for the TF
  dbExecute(con, "DROP TABLE IF EXISTS filtered_tobias;")
  dbExecute(con, paste("CREATE TEMP TABLE filtered_tobias AS
                        SELECT * FROM Tobias 
                        WHERE TF_gene_name = '", toupper(TF_name.tobias), "';", sep = ""))
  
  # Build the query dynamically
  targetQuery <- paste("
    SELECT 
      gm.gene_name,
      gm.ensg_id,
      tb.TFBS_name,
      li.feature AS feature,
      GA1_2_bound AS GA1_2_bound,
      GL1_2_bound AS GL1_2_bound,
      tb.PRO1_2_score,
      tb.CO1_2_score,
      tb.GA1_2_score,
      tb.GL1_2_score,
      tb.GA1_2_score - tb.PRO1_2_score AS PRO1_2_GA_increase,
      tb.GL1_2_score - tb.PRO1_2_score AS PRO1_2_GL_increase,
      tb.GA1_2_score / NULLIF(tb.PRO1_2_score, 0) AS PRO1_2_GA_FT_ratio,
      tb.GL1_2_score / NULLIF(tb.PRO1_2_score, 0) AS PRO1_2_GL_FT_ratio,
      tb.w_mean_cons,
      li.zscore,
      ABS(((li.end+li.start)/2) - gm.TSS_start) AS linkage_dist
    ", if (include_extra_columns) {"
      , gm.TSS_start, gm.chromosome, gm.strand, gm.exon_count
    "}, "
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
      tb.TF_gene_name = '", toupper(TF_name.tobias), "'
  ", sep = "")
  
  # Get the query result
  target.data <- as_tibble(dbGetQuery(con, targetQuery))
  
  return(target.data)
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
  
  # Keep the correlation data in pos.targets before grouping by ensg_id
  pos.targets <- pos.targets %>%
    group_by(ensg_id) %>%
    summarize(linkage_kbp = min(linkage_kbp), 
              TF.cor.PRO.GA = first(TF.cor.PRO.GA),  # Keep the correlation data for GA
              TF.cor.PRO.GL = first(TF.cor.PRO.GL))  # Keep the correlation data for GL
  
  # Create bin_category based on linkage_kbp
  pos.targets <- pos.targets %>%
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
  
  # Extract correlation data for PRO1_2 to GA and PRO1_2 to GL
  PRO1_2_GA_correlation <- pos.targets$TF.cor.PRO.GA # Correlation values for GA
  PRO1_2_GL_correlation <- pos.targets$TF.cor.PRO.GL # Correlation values for GL
  
  # Create row annotation (binary, variability + PRO1_2 GA/GL correlations)
  row_annotation_left <- rowAnnotation(
    Causal_Target = anno_simple(as.integer(rownames(data2plot) %in% gene_list), col = c("0" = "white", "1" = "green")),
    Variability = anno_simple(gene_variability, col = colorRamp2(c(min(gene_variability), max(gene_variability)), c("blue", "red")))
  )
  
  # Right side row annotation for PRO1_2 to GA and GL correlations
  row_annotation_right <- rowAnnotation(
    PRO1_2_GA_Correlation = anno_simple(PRO1_2_GA_correlation, col = colorRamp2(c(min(PRO1_2_GA_correlation), max(PRO1_2_GA_correlation)), c("blue", "red"))),
    PRO1_2_GL_Correlation = anno_simple(PRO1_2_GL_correlation, col = colorRamp2(c(min(PRO1_2_GL_correlation), max(PRO1_2_GL_correlation)), c("blue", "red")))
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
    left_annotation = row_annotation_left,        # Add left annotation to the heatmap (causal target + variability)
    right_annotation = row_annotation_right       # Add right annotation for PRO1_2 GA and GL correlations
  )
  print(h1)
  dev.off()
}

plot_heatmap_with_variability_v2 <- function(
    combined_tibble,      # Dataframe with the filtered target information from process_TF_targets_with_TE_v2
    rna_data,             # Matrix of RNA expression data
    gene_id2name,         # A named vector to map gene IDs to gene names
    target_RNA_threshold = 1.2, # RNA expression filter threshold
    output_file = "heatmap_with_variability.pdf", # Output file for the heatmap
    column_samples = c("PRO1","PRO2","CO1","CO2","GA1","GA2","GA3","GA4","GA5","GA6","GL1","GL2","GL3","GL4","GL5"), # Column names (samples) to plot
    bin_ranges = c("0-5 kbp", "> 5 kbp - 50 kbp", "> 50 kbp") # Bin ranges for splitting rows
) {
  
  # Keep the correlation data in pos.targets before grouping by ensg_id
  pos.targets <- combined_tibble %>%
    group_by(ensg_id) %>%
    reframe(linkage_kbp = min(linkage_kbp, na.rm = TRUE), 
            TF.cor.PRO.GA = first(TF.cor.PRO.GA, na.rm = TRUE),  # Keep the correlation data for GA
            TF.cor.PRO.GL = first(TF.cor.PRO.GL, na.rm = TRUE))  # Keep the correlation data for GL
  
  # Create bin_category based on linkage_kbp
  pos.targets <- pos.targets %>%
    mutate(bin_category = case_when(
      linkage_kbp <= 5 ~ bin_ranges[1],
      linkage_kbp > 5 & linkage_kbp <= 50 ~ bin_ranges[2],
      linkage_kbp > 50 ~ bin_ranges[3]
    ))
  
  # Filter expressed genes for heatmap
  exp.filt <- rownames(rna_data[rowSums(rna_data > target_RNA_threshold, na.rm = TRUE) > 0,])
  pos.targets <- pos.targets %>% filter(ensg_id %in% exp.filt)
  
  # Check if any targets remain after RNA expression filtering
  if (nrow(pos.targets) == 0) {
    stop("No targets remaining after RNA expression filtering.")
  }
  
  # Data to plot (log transformed RNA expression)
  data2plot <- log1p(rna_data[pos.targets$ensg_id, column_samples])
  
  # Set gene names as row names
  rownames(data2plot) <- sapply(rownames(data2plot), function(g) { gene_id2name[[g]] })
  
  # Calculate variability (standard deviation) per gene
  gene_variability <- apply(data2plot, 1, sd, na.rm = TRUE)
  
  # Ensure the row order of data2plot and pos.targets matches, add variability to pos.targets
  pos.targets <- pos.targets %>%
    mutate(gene_name = sapply(ensg_id, function(g) { gene_id2name[[g]] })) %>%
    filter(gene_name %in% rownames(data2plot)) %>%
    arrange(factor(gene_name, levels = rownames(data2plot))) %>%
    mutate(variability = gene_variability[match(gene_name, rownames(data2plot))])
  
  # Ensure that gene_variability matches the row names of data2plot
  gene_variability <- gene_variability[match(rownames(data2plot), names(gene_variability))]
  
  # Handle cases where all variability values are missing
  if (all(is.na(gene_variability))) {
    stop("All variability values are missing. Cannot generate heatmap.")
  }
  
  # Calculate average variability per bin_category
  avg_variability_by_bin <- pos.targets %>%
    group_by(bin_category) %>%
    summarise(avg_variability = mean(variability, na.rm = TRUE))
  
  # Append average variability to bin_category names
  pos.targets <- pos.targets %>%
    left_join(avg_variability_by_bin, by = "bin_category") %>%
    mutate(bin_category_label = paste0(bin_category, " (Avg Var: ", round(avg_variability, 3), ")"))
  
  # Extract correlation data for PRO1_2 to GA and PRO1_2 to GL
  PRO1_2_GA_correlation <- pos.targets$TF.cor.PRO.GA # Correlation values for GA
  PRO1_2_GL_correlation <- pos.targets$TF.cor.PRO.GL # Correlation values for GL
  
  # Create row annotation (binary, variability + PRO1_2 GA/GL correlations)
  row_annotation_left <- rowAnnotation(
    Variability = anno_simple(gene_variability, 
                              col = colorRamp2(c(min(gene_variability, na.rm = TRUE), 
                                                 max(gene_variability, na.rm = TRUE)), c("blue", "red")))
  )
  
  # Right side row annotation for PRO1_2 to GA and GL correlations
  row_annotation_right <- rowAnnotation(
    PRO1_2_GA_Correlation = anno_simple(PRO1_2_GA_correlation, 
                                        col = colorRamp2(c(min(PRO1_2_GA_correlation, na.rm = TRUE), 
                                                           max(PRO1_2_GA_correlation, na.rm = TRUE)), 
                                                         c("blue", "red"))),
    PRO1_2_GL_Correlation = anno_simple(PRO1_2_GL_correlation, 
                                        col = colorRamp2(c(min(PRO1_2_GL_correlation, na.rm = TRUE), 
                                                           max(PRO1_2_GL_correlation, na.rm = TRUE)), 
                                                         c("blue", "red")))
  )
  
  # Color scale for the heatmap
  color_scale <- colorRamp2(c(min(data2plot, na.rm = TRUE), max(data2plot, na.rm = TRUE)), c("white", "red"))
  
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
    left_annotation = row_annotation_left,        # Add left annotation to the heatmap (variability)
    right_annotation = row_annotation_right       # Add right annotation for PRO1_2 GA and GL correlations
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

# Helper function to handle null values
`%||%` <- function(a, b) if (!is.null(a)) a else b

process_TF_targets_with_TE_v2 <- function(input_tibble_with_te, 
                                          TF, 
                                          zscore_threshold = NULL, 
                                          pval_adj_threshold = NULL, 
                                          w_mean_cons_threshold = NULL, 
                                          pct_1_threshold = 0.5, 
                                          pct_2_threshold = 0.5, 
                                          dynamic_thresholds = TRUE, 
                                          focus_lineage = c("GA", "GL", "both"),
                                          regulation_type = c("Activation", "Repression", "Unknown"),
                                          GA_GL_bound_required = TRUE, 
                                          custom_filter = NULL, 
                                          causality = c("TF->Target", "Target->TF", "Bidirectional", "No significant TE")) {
  
  # Set dynamic thresholds based on quantiles if not provided
  if (dynamic_thresholds) {
    zscore_threshold <- zscore_threshold %||% quantile(input_tibble_with_te$zscore, 0.75, na.rm = TRUE) # Use 75th percentile
    pval_adj_threshold <- pval_adj_threshold %||% 0.05 # Use default 0.05 if not provided
    w_mean_cons_threshold <- w_mean_cons_threshold %||% quantile(input_tibble_with_te$w_mean_cons, 0.75, na.rm = TRUE) # Use 75th percentile
    pct_1_threshold <- pct_1_threshold %||% quantile(input_tibble_with_te$pct.1, 0.75, na.rm = TRUE) # Use 75th percentile
    pct_2_threshold <- pct_2_threshold %||% quantile(input_tibble_with_te$pct.2, 0.75, na.rm = TRUE) # Use 75th percentile
  }
  
  # Filter based on the dynamically calculated or provided thresholds
  filtered_targets <- input_tibble_with_te %>%
    filter(TF == !!TF) %>%
    filter(zscore > zscore_threshold &
             p_val_adj < pval_adj_threshold &
             w_mean_cons > w_mean_cons_threshold)
  
  # Debug: Print how many rows are left after this filter
  cat("Number of target genes after zscore, p_val_adj, and w_mean_cons filtering:", nrow(filtered_targets), "\n")
  
  # Apply lineage-specific filtering if required
  if (focus_lineage == "GA") {
    filtered_targets <- filtered_targets %>%
      filter(GA1_2_bound == 1)
  } else if (focus_lineage == "GL") {
    filtered_targets <- filtered_targets %>%
      filter(GL1_2_bound == 1)
  }
  
  # Debug: Print how many rows are left after lineage filter
  cat("Number of target genes after lineage filtering:", nrow(filtered_targets), "\n")
  
  # Apply custom filter if provided
  if (!is.null(custom_filter)) {
    filtered_targets <- filtered_targets %>%
      filter(eval(parse(text = custom_filter)))
  }
  
  # Debug: Print how many rows are left after custom filter
  cat("Number of target genes after custom filtering:", nrow(filtered_targets), "\n")
  
  # Filter based on regulation type
  filtered_targets <- filtered_targets %>%
    mutate(regulation = case_when(
      zscore > 0 & avg_log2FC > 0.1 ~ "Activation",  # TF activates the gene
      zscore < 0 & avg_log2FC < -0.1 ~ "Repression", # TF represses the gene
      TRUE ~ "Unknown"                               # Mixed or weak signal
    )) %>%
    filter(regulation %in% regulation_type)
  
  # Debug: Print how many rows are left after regulation filtering
  cat("Number of target genes after regulation type filtering:", nrow(filtered_targets), "\n")
  
  # Filter based on transfer entropy results and causality
  filtered_targets <- filtered_targets %>%
    mutate(
      causality_present_X_Y = ifelse(X_Y_te_p_value < 0.05 & X_Y_te_coef > 0.02, TRUE, FALSE), # X->Y TE significant and strong
      causality_present_Y_X = ifelse(Y_X_te_p_value < 0.05 & Y_X_te_coef > 0.02, TRUE, FALSE), # Y->X TE significant and strong
      TE_relationship = case_when(
        causality_present_X_Y & !causality_present_Y_X ~ "TF->Target",
        causality_present_Y_X & !causality_present_X_Y ~ "Target->TF",
        causality_present_X_Y & causality_present_Y_X ~ "Bidirectional",
        TRUE ~ "No significant TE"
      )
    ) %>%
    filter(TE_relationship %in% causality)
  
  # Debug: Print how many rows are left after causality filtering
  cat("Number of target genes after causality filtering:", nrow(filtered_targets), "\n")
  
  # Filter based on pct.1 and pct.2 thresholds for GA and GL expression
  if (focus_lineage == "GA") {
    filtered_targets <- filtered_targets %>%
      filter(pct.1 > pct_1_threshold)
  } else if (focus_lineage == "GL") {
    filtered_targets <- filtered_targets %>%
      filter(pct.2 > pct_2_threshold)
  }
  
  # Debug: Print how many rows are left after pct filtering
  cat("Number of target genes after pct filtering:", nrow(filtered_targets), "\n")
  
  # Keep GA/GL bound condition if required
  if (GA_GL_bound_required) {
    filtered_targets <- filtered_targets %>%
      filter(GA1_2_bound == 1 | GL1_2_bound == 1)
  }
  
  # Debug: Final number of rows
  cat("Final number of target genes after all filters:", nrow(filtered_targets), "\n")
  
  return(filtered_targets)
}

# Function to create barplot showing percentage of target genes in expression categories
plot_expression_categories <- function(
    combined_tibble,  # Dataframe with target gene information
    rna_data,         # RNA expression data with genes as rows and samples as columns
    cell_groups = c("GA", "GL")  # Cell groups to compare (GA and GL)
) {
  
  # Ensure target genes (ensg_id) are present as rows in the RNA data
  valid_targets <- combined_tibble %>%
    filter(ensg_id %in% rownames(rna_data))
  
  if (nrow(valid_targets) == 0) {
    stop("No valid target genes found in rna_data.")
  }
  
  # Extract expression data for the valid target genes
  expression_data <- rna_data[valid_targets$ensg_id, , drop = FALSE]
  
  # Reshape the expression data for categorization, grouping by GA and GL cell groups
  plot_data <- expression_data %>%
    as.data.frame() %>%
    mutate(ensg_id = rownames(expression_data)) %>%
    pivot_longer(cols = -ensg_id, names_to = "Sample", values_to = "Expression") %>%
    left_join(valid_targets, by = "ensg_id") %>%
    mutate(
      group = case_when(
        grepl("^GA", Sample) ~ "GA Group",
        grepl("^GL", Sample) ~ "GL Group"
      )
    ) %>%
    filter(!is.na(group))  # Remove any rows without GA or GL group
  
  # Categorize expression values into three bins: <0.5, 0.5-1.5, >1.5
  plot_data <- plot_data %>%
    mutate(expression_category = case_when(
      Expression < 0.5 ~ "< 0.5",
      Expression >= 0.5 & Expression <= 1.5 ~ "0.5 - 1.5",
      Expression > 1.5 ~ "> 1.5"
    ))
  
  # Calculate the percentage of genes in each category for GA and GL
  percentage_data <- plot_data %>%
    group_by(group, expression_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(percentage = (count / sum(count)) * 100)
  
  # Create the barplot
  ggplot(percentage_data, aes(x = expression_category, y = percentage, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      x = "Expression Category",
      y = "Percentage of Target Genes (%)",
      title = "Percentage of Target Genes in Expression Categories (GA vs GL)"
    ) +
    theme_minimal()
}


# Function to create boxplot of target gene expression levels in GA and GL cell groups
plot_expression_by_cell_group <- function(
    combined_tibble,  # Dataframe with target gene information
    rna_data,         # RNA expression data with genes as rows and samples as columns
    cell_groups = c("GA", "GL")  # Cell groups to compare (GA and GL)
) {
  
  # Ensure target genes (ensg_id) are present as rows in the RNA data
  valid_targets <- combined_tibble %>%
    filter(ensg_id %in% rownames(rna_data))
  
  if (nrow(valid_targets) == 0) {
    stop("No valid target genes found in rna_data.")
  }
  
  # Extract expression data for the valid target genes
  expression_data <- rna_data[valid_targets$ensg_id, , drop = FALSE]
  
  # Reshape the expression data for plotting, grouping by GA and GL cell groups
  plot_data <- expression_data %>%
    as.data.frame() %>%
    mutate(ensg_id = rownames(expression_data)) %>%
    pivot_longer(cols = -ensg_id, names_to = "Sample", values_to = "Expression") %>%
    left_join(valid_targets, by = "ensg_id") %>%
    mutate(
      group = case_when(
        grepl("^GA", Sample) ~ "GA Group",
        grepl("^GL", Sample) ~ "GL Group"
      )
    ) %>%
    filter(!is.na(group))  # Remove any rows without GA or GL group
  
  # Create the boxplot
  ggplot(plot_data, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot(fill = "grey") +
    labs(
      x = "Cell Group",
      y = "Gene Expression Level",
      title = "Target Expression in GA and GL lineages"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18,angle = 45, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 16)) + xlab("")
}

# Function to plot the barplot of zscore counts
plot_zscore_bar <- function(data) {
  # Create a summary table with counts for zscore > 0 and zscore < 0
  zscore_summary <- data %>%
    summarize(
      greater_than_zero = sum(zscore > 0),
      less_than_zero = sum(zscore < 0)
    ) %>%
    pivot_longer(cols = everything(), names_to = "zscore_category", values_to = "count")
  
  # Create the barplot using ggplot2
  ggplot(zscore_summary, aes(x = zscore_category, y = count, fill = zscore_category)) +
    geom_bar(stat = "identity",color = "grey", fill = "grey") +
    scale_x_discrete(labels = c("greater_than_zero" = "Pos. link", "less_than_zero" = "Neg. link")) +
    labs(
      x = "Link direction",
      y = "Count",
      title = ""
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18,angle = 45, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 16)) + xlab("")
}

# Function to plot the barplot of linkage_kbp grouped into categories with closest linkage counted only once per gene
plot_linkage_kbp_bar <- function(data, count_closest_only = TRUE) {
  
  # Filter out rows with NA in linkage_kbp
  data <- data %>%
    filter(!is.na(linkage_kbp))
  
  # If count_closest_only is TRUE, assign each gene to the category of its closest linkage
  if (count_closest_only) {
    data <- data %>%
      group_by(ensg_id) %>%
      slice_min(linkage_kbp, n = 1, with_ties = FALSE) %>%  # Keep only the closest linkage, without ties
      ungroup()
  }
  
  # Create categories based on linkage_kbp ranges
  linkage_summary <- data %>%
    mutate(linkage_category = case_when(
      linkage_kbp <= 5 ~ "0-5 kbp",
      linkage_kbp > 5 & linkage_kbp <= 50 ~ "5-50 kbp",
      linkage_kbp > 50 ~ "> 50 kbp"
    )) %>%
    # Group by the categories and count how many fall into each category
    group_by(linkage_category) %>%
    summarize(count = n()) %>%
    # Reorder the linkage categories for proper ordering on the x-axis
    mutate(linkage_category = factor(linkage_category, levels = c("0-5 kbp", "5-50 kbp", "> 50 kbp")))
  
  # Create the barplot using ggplot2
  ggplot(linkage_summary, aes(x = linkage_category, y = count, fill = linkage_category)) +
    geom_bar(stat = "identity", color = "grey", fill = "grey") +
    labs(
      x = "Linkage Distance Category (kbp)",
      y = "Count",
      title = ""
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 16)) + 
    xlab("")
}
# Simple function to plot boxplot of gene expression variability by linkage distance category
plot_simple_variability_boxplot <- function(
    combined_tibble, 
    rna_data, 
    bin_ranges = c("0-5 kbp", "5-50 kbp", "> 50 kbp"),
    add_p_value = FALSE  # New argument to toggle p-value testing
) {
  
  # Create categories based on linkage_kbp ranges
  pos.targets <- combined_tibble %>%
    filter(!is.na(linkage_kbp)) %>%  # Remove rows with NA values in linkage_kbp
    mutate(bin_category = case_when(
      linkage_kbp <= 5 ~ bin_ranges[1],
      linkage_kbp > 5 & linkage_kbp <= 50 ~ bin_ranges[2],
      linkage_kbp > 50 ~ bin_ranges[3]
    )) %>%
    filter(ensg_id %in% rownames(rna_data)) %>%
    ungroup()
  
  # Fetch and log-transform RNA expression data
  expression_data <- log1p(rna_data[pos.targets$ensg_id, ])
  
  # Calculate standard deviation (variability) per gene
  gene_variability <- apply(expression_data, 1, sd)
  
  # Add gene variability to the pos.targets tibble
  pos.targets <- pos.targets %>%
    mutate(variability = gene_variability[match(ensg_id, rownames(expression_data))])
  
  # Ensure correct x-axis order by defining factor levels
  pos.targets <- pos.targets %>%
    mutate(bin_category = factor(bin_category, levels = bin_ranges))
  
  # Define comparisons (what pairs to test against each other)
  comparisons <- list(
    c("0-5 kbp", "5-50 kbp"),
    c("0-5 kbp", "> 50 kbp"),
    c("5-50 kbp", "> 50 kbp")
  )
  
  # Start building the boxplot
  p <- ggplot(pos.targets, aes(x = bin_category, y = variability, fill = bin_category)) +
    geom_boxplot(fill = "grey") +
    labs(
      x = "Linkage Distance Category",
      y = "Expression Variability (Standard Deviation)",
      title = ""
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 18,angle = 45, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 16))  + xlab("")
  
  # Add p-value comparisons if add_p_value is TRUE
  if (add_p_value) {
    p <- p + stat_compare_means(comparisons = comparisons, method = "t.test", 
                                label = "p.format", tip.length = 0.03)  # Add p-values with connecting lines
  }
  
  # Print the plot
  print(p)
}

# Function to create scatter plot of Tal1 vs target gene expression
plot_tal1_vs_target_scatter <- function(
    tal1_ensg_id,        # The ensg_id of Tal1 (e.g., "ENSMUSG00000012356")
    target_genes,        # A vector of target gene ensg_ids
    rna_data,            # RNA expression matrix
    columns = c("CO2", "GA1")  # Columns to use for plotting
) {
  
  # Ensure Tal1 is present in the RNA data
  if (!tal1_ensg_id %in% rownames(rna_data)) {
    stop("Tal1 ensg_id not found in rna_data.")
  }
  
  # Ensure the target genes are present in the RNA data
  valid_targets <- target_genes[target_genes %in% rownames(rna_data)]
  
  if (length(valid_targets) == 0) {
    stop("No valid target genes found in rna_data.")
  }
  
  # Create an empty data frame to store results
  plot_data <- data.frame()
  
  # Extract Tal1 expression values for the selected columns
  tal1_expression <- rna_data[tal1_ensg_id, columns]
  
  # Loop over valid target genes and create scatter plot data
  for (target_gene in valid_targets) {
    target_expression <- rna_data[target_gene, columns]
    
    # Combine Tal1 and target gene expression into one data frame
    temp_data <- data.frame(
      Tal1_expression = tal1_expression,
      Target_expression = target_expression,
      Condition = columns,
      Target_gene = target_gene
    )
    
    # Append to the final plot data
    plot_data <- rbind(plot_data, temp_data)
  }
  
  # Create the scatter plot
  ggplot(plot_data, aes(x = Tal1_expression, y = Target_expression, color = Condition)) +
    geom_point() +
    facet_wrap(~ Target_gene, scales = "free_y") +  # Create separate panels for each target gene
    labs(
      x = "Tal1 Expression (CO2 / GA1)",
      y = "Target Gene Expression (CO2 / GA1)",
      title = "Scatter Plot of Tal1 vs Target Gene Expression"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 18,angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 16)) + xlab("")
}

# Function to create boxplot of target gene expression levels by linkage distance category
plot_expression_by_distance_category <- function(
    combined_tibble,  # Dataframe with target gene information
    rna_data,         # RNA expression data with genes as columns and samples as rows
    bin_ranges = c("0-5 kbp", "5-50 kbp", "> 50 kbp")  # Distance bins
) {
  
  # Filter combined_tibble for valid target genes that are present as columns in the RNA data
  pos.targets <- combined_tibble %>%
    filter(ensg_id %in% colnames(rna_data) & !is.na(linkage_kbp)) %>%
    group_by(ensg_id) %>%
    slice(1) %>%  # Keep only the first occurrence of each ensg_id to avoid duplicates
    ungroup() %>%
    mutate(bin_category = case_when(
      linkage_kbp <= 5 ~ bin_ranges[1],
      linkage_kbp > 5 & linkage_kbp <= 50 ~ bin_ranges[2],
      linkage_kbp > 50 ~ bin_ranges[3]
    ))
  
  # Extract the expression data for the filtered target genes
  expression_data <- rna_data[, pos.targets$ensg_id, drop = FALSE]
  
  # Reshape the data for plotting, excluding the Sample column
  plot_data <- expression_data %>%
    as.data.frame() %>%
    mutate(Sample = rownames(rna_data)) %>%
    pivot_longer(cols = -Sample, names_to = "ensg_id", values_to = "Expression") %>%
    left_join(pos.targets, by = "ensg_id")
  
  # Ensure correct x-axis order by setting factor levels for bin_category
  plot_data <- plot_data %>%
    mutate(bin_category = factor(bin_category, levels = bin_ranges))
  
  # Create the boxplot
  ggplot(plot_data, aes(x = bin_category, y = Expression, fill = bin_category)) +
    geom_boxplot() +
    labs(
      x = "Linkage Distance Category",
      y = "Gene Expression Level",
      title = "Target Gene Expression by Linkage Distance"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}