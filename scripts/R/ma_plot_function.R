##################
## Code for MA plots
######################



# Run: plot_ma_by_rep(condition_combis = generate_condition_combinations(dds, "ctrl"))

# Takes in dds and condition of interest. Generates all pairwise comparisons possible, returns it in a list. 
# Use this as input for plot_ma_by_rep

generate_condition_combinations <- function(dds, condition_of_interest){
  
  # Subset samples
  samples = colData(dds)[colData(dds)$condition == condition_of_interest, ]$sample_id
  
  # Generate combinations
  condition_combis <- combn(samples, m = 2, simplify = FALSE)
  
  return(condition_combis)
}

plot_ma_by_rep <- function(condition_combis, ncol = 2){
  
  ma_plots <- lapply(condition_combis, function(comparison){
    
    x = counts(dds, normalized = TRUE)[, as.character(comparison[1])]
    y = counts(dds, normalized = TRUE)[, as.character(comparison[2])]
    
    M = log2(x) - log2(y)
    A = (log2(x) + log2(y)) / 2
    df = data.frame(gene_id = names(x), M = M, A = A)
    
    # Plot graph  
    suppressWarnings({
      
      alpha <- ifelse(abs(df$M) < 1.5, 0.02, 0.2)
      p1 <- ggplot(df, aes(x = A, y = M)) + 
        # Draw stuff
        geom_point(size = 1.5, alpha = 0.2
                   #alpha = alpha
        ) + 
        stat_smooth(se = FALSE, method = "loess", color = "red3") + 
        #gghighlight::gghighlight(abs(M) > 1.5) +
        geom_hline(yintercept = 0, 
                   colour = "blue3", 
                   linetype = "dashed") + 
        geom_hline(yintercept = 1.5, color = "darkgreen", linetype = "dashed") +
        geom_hline(yintercept = -1.5, color = "darkgreen", linetype = "dashed") +
        
        # Scales, titles, etc
        scale_y_continuous(limits = c(-6, 6),
                           breaks = c(-6, -1.5, 0, 1.5, 6)) +
        labs(title = paste(comparison[1], "vs", comparison[2], sep = " "),
             subtitle = paste("Median abs diff =", 
                              round(median(abs(x - y)), digits = 2), sep = " ")
             ) + 
        scale_alpha_continuous(guide = FALSE) +
        theme_classic()
      
    })
    
  })
  
  # Plot
  suppressWarnings(
    suppressMessages(
      # if only 1 plot, no need columns
      if (length(ma_plots) == 1){
        print(ma_plots)
      } else {
        gridExtra::grid.arrange(
          grobs = ma_plots, ncol = ncol,
        )
      }
    )
  )
  
}