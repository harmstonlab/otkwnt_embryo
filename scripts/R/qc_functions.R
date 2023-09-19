# QC functions

make_pca <- function(rld, intgroup,  
                     title = "title",         
                     xlimits = c(-30, 30), 
                     ylimits = c(-15, 15),
                     manual_colors = NA,
                     label = FALSE, 
                     rearrange_intgroup = FALSE,
                     intgroup_order){
  # Calculations
  ntop = 500
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  data <- plotPCA(rld, intgroup = intgroup, 
                  returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"), digits = 2)
  
  if(rearrange_intgroup == TRUE){
    data[[intgroup]] <- factor(data[[intgroup]], 
                               levels = intgroup_order)
  }

  
  # Plot 
  pca_plot <- ggplot(as.data.frame(data), aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = as.data.frame(data)[[intgroup]]), 
               size = 2, alpha = 0.8) +
    labs(title = title,
         subtitle = paste0("By ", intgroup), 
         colour = intgroup) +
    
    # Add scale annotations
    scale_x_continuous(paste0("PC1: ",percentVar[1],"% variance"), 
                       limits = xlimits) +
    scale_y_continuous(paste0("PC2: ",percentVar[2],"% variance"), 
                       limits = ylimits) +
    # Default: paired
    # scale_color_brewer(palette = "Paired") +
    # specify custom colors
    scale_color_manual(name = intgroup, 
                       values = manual_colors) +
    
    # Make 1 unit on x-axis equal to 1 unit on y-axis
    coord_fixed(ratio = 1) +
    theme_classic()
  
  if(label == TRUE){
    pca_plot <- pca_plot + 
      geom_text_repel(data = data, aes(PC1,PC2, label = name), 
                      hjust = 0.5, box.padding = 0.5, size = 3,
                      max.overlaps = Inf)
    return(pca_plot)
  } else {
    return(pca_plot)
  }
  
}
