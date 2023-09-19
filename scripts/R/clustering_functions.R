



cluster_kmeans <- function(rld_z, nclust, plot_sil = FALSE){
  
  ## Perform k means clustering
  set.seed(1)
  nclust = nclust
  kmeans_coef = kmeans(rld_z, nclust, nstart = 1000, iter.max = 50)
  
  ## Add silhouette plots for diagnosis 
  if(plot_sil == TRUE){
    # Set colors
    c8 <- c("#EB6424","#FA9500", "#FFBB00",
            "#FFEE00", "#CCBA00",
            "#B99900", "#AA9900", "#AAAA99") 
    
    # Calculate distance
    d <- dist(rld_z)
    
    # Plot
    plot(cluster::silhouette(kmeans_coef$cluster, d),
         col= c8[1:nclust],
         border = NA, 
         main = paste("Silhouette plot, k = ", nclust)
    )
    
    return(kmeans_coef)
    
  } else {
    return(kmeans_coef)
  }
}



plot_kmeans_heatmap <- function(rld_z, k_coef,
                                sample_order, ...){
  
  # Arrange samples in correct order
  #order_samples <- as.factor(sample_order)
  
  #rld_z <- rld_z[ ,order_samples]
  
  # Set up heatmap
  # Thresholding it to 3 standard deviations away from the median. 
  thr = 3
  rld_z[rld_z > thr] = thr
  rld_z[rld_z < -thr] = -thr
  
  # Sort out color scheme
  paletteLength = 20 
  breaksList <- c(seq(-thr, 0, length.out = ceiling(paletteLength/2) + 1), 
                  seq(thr/paletteLength, thr, length.out = floor(paletteLength/2)))
  
  color = c(colorRampPalette(c("mediumblue", "white"))(14), 
            colorRampPalette(c("white", "firebrick2"))(14))
  breaksList = seq(-3, 3, length.out = 29)
  
  cs = k_coef$cluster
  #cs <- factor(cs, levels = c(2, 1))
  
  # order things in correct order
  z.toplot = rld_z[order(cs), sample_order]
  
  
  # Heatmap according to clusters
  heatmap_cl = pheatmap::pheatmap(z.toplot, 
                                  breaks = breaksList, 
                                  cluster_col = FALSE,
                                  cluster_rows = FALSE, 
                                  show_rownames = FALSE,
                                  show_colnames = TRUE, 
                                  #color = color,
                                  #annotation_col = annotation,
                                  #annotation_colours = anno_colours,
                                  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(29),
                                  fontsize_col = 12, 
                                  legend = TRUE,
                                  border_color = NA, 
                                  main = paste("Heatmap"),
                                  angle_col = 45,
                                  ...
  )
  
  heatmap_cl
  
}


get_cluster_genes <- function(k_coef, nclust){
  
  # Loop to obtain genes for each cluster
  kmeans_cl = c()
  
  for (i in 1:nclust){
    kmeans_cl[[i]] = k_coef$cluster[
      which(k_coef$cluster == i)
    ]
  }
  
  return(kmeans_cl)
}


# Run: zscore_boxplot(kmeans_cl_k6, clust_num = 1, sample_order)

zscore_boxplot <- function(kmeans_cl, clust_num, 
                           sample_order = sample_order){
  # Get the genes
  c3_genes <- names((kmeans_cl)[[clust_num]])
  print(paste0(length(c3_genes), " genes in cluster ", clust_num))
  
  # Plot the normalized counts
  c3_all <- rld_z[c3_genes, ] %>% 
    melt()
  colnames(c3_all) <- c("gene_id", "sample", "vst_z")
  # Reorder samples
  sample_order
  c3_all$sample <- factor(c3_all$sample, levels = sample_order)
  
  # Plot the zscores
  zscore_plot <- ggplot(c3_all, aes(x = sample, y = vst_z)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
    geom_boxplot() +
    scale_y_continuous(limits = c(-4, 4)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste0("Zscore boxplot, cluster ", clust_num),
         x = "", 
         y = "z-score") 
  
  return(zscore_plot)
  
}


# Run: zscore_boxplot(kmeans_cl_k6, clust_num = 1, sample_order)

zscore_boxcondition <- function(kmeans_cl, clust_num, 
                           condition_order = condition_order){
  # Get the genes
  c3_genes <- names((kmeans_cl)[[clust_num]])
  print(paste0(length(c3_genes), " genes in cluster ", clust_num))
  
  # Plot the normalized counts
  c3_all <- rld_z[c3_genes, ] %>% 
    melt()
  colnames(c3_all) <- c("gene_id", "sample", "vst_z")
  # Reorder conditions
  c3_all$condition <- c3_all$sample
  c3_all$condition <- gsub("_[1-3]", "", c3_all$condition)
  c3_all$condition <- factor(c3_all$condition, levels = condition_order)
  
  # Plot the zscores
  zscore_plot <- ggplot(c3_all, aes(x = condition, y = vst_z)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
    geom_boxplot() +
    scale_y_continuous(limits = c(-4, 4)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste0("Zscore boxplot, cluster ", clust_num),
         x = "", 
         y = "z-score") 
  
  return(zscore_plot)
  
}



### ------ clusterHeatmap ------ ###

clusterHeatmap = function(rld_z,
                          kmeans_cl, 
                          clust_num, 
                          sample_order,
                          cluster_rows = TRUE,
                          cluster_columns = FALSE, 
                          show_row_names = FALSE,
                          show_column_names = TRUE,
                          font_size = 14
) {
  
  # Set up heatmap
  # Thresholding it to 3 standard deviations away from the median. 
  thr = 3
  rld_z[rld_z > thr] = thr
  rld_z[rld_z < -thr] = -thr
  
  # Subset rld   
  rld_z_cl = rld_z[which(rownames(rld_z) %in% 
                           names(kmeans_cl[[ clust_num ]])), ]
  
  # Create a heatmap
  
  heatmap_cl = ComplexHeatmap::Heatmap(
    matrix = rld_z_cl, 
    col = colorRamp2(c(-3, 0, 3), 
                     c("mediumblue", "white", 
                       "firebrick2")
    ),
    column_title = paste0("Heatmap for cluster ", clust_num),
    column_title_side = "top",
    name = "z-score",
    ## Top annotation (couldn't get this to display properly)
    #top_annotation = ha,
    
    ## -- Cluster rows --- ## 
    cluster_rows = cluster_rows,
    clustering_distance_rows = "euclidean", # clust by genes 
    
    ## Rows: graphic parameters  
    
    show_row_names = show_row_names,
    row_names_gp =  gpar(fontsize = 11),
    row_names_side = "right",
    row_dend_gp = gpar(fontsize = 12),  
    row_dend_width = unit(20, "mm"),
    
    ## -- Cluster columns -- ##
    cluster_columns = cluster_columns,
    
    ## Columns: graphic parameters
    show_column_names = show_column_names, 
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 14), 
    column_names_side = "bottom",
    column_dend_gp = gpar(fontsize = 12),
    column_dend_side = "top",
    column_dend_height = unit(10, "mm"),
    
    column_order = factor(sample_order, levels = sample_order)
  )
  heatmap_cl
  
}






