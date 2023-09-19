#################

# DE-functions

##################

# Description: Functions used in 02_DE for a typical differential expression analysis. 
# Author: Qian Hui TAN
# Date last modified: 18-Apr-2023



## This function is to write out the files of differentially expressed genes
## It takes in a DESeqResults object, and two strings -- the numerator and denominator used in the analysis -- and writes out csv files

write_files <- function(results, numerator, denominator, output_directory = output_dir){
  # these are all the genes that are differentially expressed between the two conditions, not just the significant ones
  write.csv(results, paste0(output_directory,numerator,"_",denominator,"_all.csv"), row.names = TRUE, col.names = TRUE, sep = " ")
  
  # these are the genes that are significantly differentially expressed by FDR 10% and abs(log2fc) > log2(1.5)
  sig_padj_genes <- results[!is.na(results$padj),]
  sig_padj_genes <- sig_padj_genes[sig_padj_genes$padj < 0.1,]
  sig_padj_fc_genes <- sig_padj_genes[abs(sig_padj_genes$log2FoldChange) > lfc.threshold,]
  write.csv(sig_padj_fc_genes, paste0(output_dir,numerator,"_",denominator,"_significant.csv"), row.names = TRUE, col.names = TRUE)
}


## This function generates the MA plots with significant changes above the threshold coloured in red and significant changes below the threshold coloured in blue
## It takes in a DESeqResults object, uses the plotMA function from DESeq2 to obtain the necessary data frame to plot

generate_ma <- function(results){
  df <- DESeq2::plotMA(results, ylim = c(-10,10), colSig = "red", returnData = TRUE)
  plot <- df %>%
    mutate(signif = ifelse(lfc > lfc.threshold & isDE == TRUE, "U", 
                           ifelse(lfc < -lfc.threshold & isDE == TRUE, "D", "N"))) %>%
    ggplot(aes(x=mean, y=lfc, colour = signif)) + 
    geom_point(size = 1.5, alpha = 0.8) + 
    theme_classic() + 
    geom_hline(yintercept=0, colour="grey40", lwd = 1) + 
    #stat_smooth(se = FALSE, method = "loess", color = "red3") + 
    theme_classic() + 
    scale_colour_manual(values=c("#4575b4","#a3a3a3","#d73027"), labels = c("D", "N", "U")) +
    ylim(c(-10,10)) +
    theme(legend.position = "none") +
    ylab("Log fold change") +
    xlab("Mean of normalized counts") +
    scale_x_log10()
  return(plot)
}

## This function plots the volcano plot
## It takes in a data frame and two strings which are used for the title of the plot

generate_volcano <- function(data_frame, numerator, denominator){
  lfc.threshold = log2(1.5)
  tmp = as.data.frame(data_frame)
  tmp$signif = ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.01, "U1", 
                      ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.05, "U2",
                             ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.1, "U3",
                                    ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.01, "D1", 
                                           ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.05, "D2",
                                                  ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.1, "D3",                                                  "N"))))))
  tmp$signif = factor(tmp$signif, c("N", "U1", "U2", "U3", "D3", "D2", "D1"))
  
  x = ggplot(data=tmp, aes(x=log2FoldChange, y=-log10(padj), colour= signif)) + geom_point(alpha=1.0, size=2.00) +
    ggtitle(paste("Volcano Plot:", numerator, "vs.", denominator)) + scale_x_continuous("log2(fold change)", limits=c(-10, 10)) +    
    scale_y_continuous("-log10(FDR)") + geom_vline(xintercept = lfc.threshold, linetype="dotdash") +
    geom_vline(xintercept = -1*(lfc.threshold), linetype="dotdash") +
    geom_hline(yintercept = -log10(0.1), colour="gray40", linetype="dotdash") +   
    geom_hline(yintercept = -log10(0.05), colour="gray40", linetype="dotdash") + 
    geom_hline(yintercept = -log10(0.01), colour="gray40", linetype="dotdash") + 
    scale_colour_manual("", values=c("#666666", "#d73027", "#f46d43", "#fdae61", "#abd9e9", "#74add1", "#4575b4" ), labels = c("N", "U1", "U2", "U3", "D3", "D2", "D1")) + theme_classic() + theme(legend.position = "none", plot.title = element_text(size = 20), axis.title=element_text(size=16,face="bold"))
  return(x)
  print(table(tmp$signif))
}




## This function gets a list of lfcShrink fold changes from a wald dds
## It takes in an nbinomWaldTest object, a character string contrast, and an ensembl.genes object containing gene biotype for annotations. 

get_dds_res <- function(wald_dds, contrast, ensembl.genes, 
                        lfcshrinktype = "apeglm", 
                        shrink = TRUE, parallel = TRUE) {
  
  # Run: get_dds_res(wald_dds, 
  #                  contrast = c("condition", "mdd", "ctrl"),
  #                  ensembl.genes = ensembl.genes,
  #                  shrink = TRUE)
  
  # Get the results
  results_c1_control = results(wald_dds, 
                               contrast = contrast,  
                               filter = rowMeans(counts(wald_dds, 
                                                        normalized = TRUE)), 
                               test = "Wald", alpha = 0.1, 
                               independentFiltering = TRUE)
  
  # Make the condition_mdd_vs_ctrl string
  deseq_coef <- paste(contrast[1], contrast[2], "vs", 
                      contrast[3], sep = "_")
  print(deseq_coef)
  
  # Shrink
  if (shrink == TRUE){
    results_c1_control <- lfcShrink(wald_dds, 
                                    coef = deseq_coef, 
                                    res = results_c1_control, 
                                    type = lfcshrinktype, parallel = TRUE)
  }
  
  # Add gene annotations
  results_c1_control$gene_biotype = ensembl.genes$gene_biotype[match(row.names(results_c1_control), ensembl.genes$gene_id)]
  results_c1_control$external_gene_name = ensembl.genes$external_gene_name[match(row.names(results_c1_control), ensembl.genes$gene_id)]
  
  return(results_c1_control)
}