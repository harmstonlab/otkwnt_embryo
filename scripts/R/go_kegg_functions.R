#################

# GO-KEGG-functions

##################

# Description: Functions used for GO and KEGG enrichments
# Author: Qian Hui TAN
# Date last modified: 19-Sept-2023


### --- GO functions --- ###

plotEGO_dm = function(clust_target_genes,
                      universe = universe, 
                      ont = "BP", 
                      title = "title", 
                      subtitle = NULL, 
                      font_size = 14){
  
  # Run KEGG enrichment
  message("Running GO for organism = drosophila melanogaster")
  
  cl_target_ego = enrichGO(gene = clust_target_genes, 
                           universe = universe,
                           OrgDb = org.Dm.eg.db,
                           keyType = 'ENSEMBL', 
                           ont = ont, 
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.1,
                           qvalueCutoff = 0.1,
                           readable = TRUE)
  
  # If no GO terms found, return warning message and a tibble with NA
  if(nrow(data.frame(cl_target_ego)) == 0) {
    warning(paste0("No GO enrichment found. Returning a NA tibble."
    ))
    return(tibble(`Description` = "NA"))
  } else {
    print(dotplot(cl_target_ego,
                  title = title, 
                  font.size = font_size))
    # Print number of enrichments found
    print(paste0(nrow(cl_target_ego), " enrichments found"))
    
    return(as_tibble(cl_target_ego))
  }
  
}


custom_ego <- function(ego_tibble,
                       interesting_pathways, 
                       title = "title", 
                       font_size = 14) {
  
  # Get relevant columns
  c2_interesting_egos <- ego_tibble%>% 
    dplyr::select("Description", "GeneRatio", "p.adjust", "Count") %>% 
    filter(Description %in% interesting_pathways) %>% 
    # Convert GeneRatio from fraction to decimal
    mutate(GeneRatio = DOSE::parse_ratio(GeneRatio)) %>% 
    arrange(GeneRatio) %>% 
    # Preserve sorted order
    mutate(Description = factor(Description, levels = unique(Description)))
  
  
  # Plot
  df_plot <- ggplot(c2_interesting_egos, aes(x = GeneRatio, y = Description,
                                             color = p.adjust)) +
    geom_point(aes(size = Count)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = title,
      y = "") +
    theme_light() +
    theme(axis.text.y = element_text(size = font_size))
  
  return(df_plot)
  
}


### --- KEGG functions --- ###


get_entrez <-function(ensembl_ids, ensembl_genes){
  entrez <- na.omit(ensembl.genes[ensembl.genes$gene_id %in% ensembl_ids, ]$entrezgene_id)
  return(entrez)
}

# Takes in an ego tibble, subsets interesting pathways. 
custom_ego_table <- function(ego_tibble, interesting_pathways) {
  ego_tbl <- ego_tibble[ego_tibble$Description %in% interesting_pathways, 
                        colnames(ego_tibble) %in% c("Description", "GeneRatio", "p.adjust", "geneID")]
  return(ego_tbl)
}

# Takes in an kegg tibble, subsets interesting pathways. 
custom_kegg_table <- function(ego_tibble, interesting_pathways) {
  
  ego_tbl <- ego_tibble
  ego_tbl$Description = gsub(" - Drosophila melanogaster (fruit fly)", "", ego_tbl$Description,
                             fixed = TRUE)
  ego_tbl <- ego_tbl[ego_tbl$Description %in% interesting_pathways, 
                     colnames(ego_tbl) %in% c("Description", "GeneRatio", "p.adjust", "geneID")]
  return(ego_tbl)
}



### This function ONLY WORKS FOR HUMAN GENES

# Run: 
## Get the entrez IDs for cluster 2
# c1_entrez <- na.omit(ensembl.genes[ensembl.genes$gene_id %in% names((kmeans_cl_k2)[[1]]), ]$entrezgene_id)
# k2_c1_kegg <- plotKEGG(c1_entrez, title = "KEGG, cluster 1")

plotKEGG_dm <- function(target_genes_entrez, title = "title", font_size = 14){
  
  # Run KEGG enrichment
  message("Running KEGG for organism = drosophila melanogaster")
  ekegg = clusterProfiler::enrichKEGG(target_genes_entrez,
                                      organism = "dme",
                                      keyType = "ncbi-geneid")
  
  # If no GO terms found, return warning message and a tibble with NA
  if(nrow(data.frame(ekegg)) == 0) {
    warning(paste0("No KEGG enrichment found. Returning a NA tibble."
    ))
    return(tibble(`Description` = "NA"))
  } else {
    
    # Print the kegg dotplot
    print(dotplot(ekegg, title = title, 
                  font.size = font_size))
    
    # Convert entrez to gene symbols
    ekegg = DOSE::setReadable(ekegg, OrgDb = "org.Dm.eg.db", keyType = "ENTREZID")
    
    # Print number of enrichments found
    print(paste0(nrow(ekegg), " enrichments found"))
    
    # Return the result as a data frame
    return(as.data.frame(ekegg))
  }
}




custom_ekegg <- function(ekegg_tibble, 
                         interesting_pathways, 
                         title = "title", 
                         font_size = 14) {
  # Remove trailing drosophila from Description
  ekegg_tibble$Description = gsub(" - Drosophila melanogaster (fruit fly)", "", 
                            fixed = TRUE, ekegg_tibble$Description)
  
  # Get relevant columns
  interesting_ekeggs  <- ekegg_tibble %>% 
    dplyr::select("Description", "GeneRatio", "p.adjust", "Count") %>% 
    filter(Description %in% interesting_pathways) %>%
    # Convert GeneRatio from fraction to decimal 
    mutate(GeneRatio = DOSE::parse_ratio(GeneRatio)) %>% 
    arrange(GeneRatio) %>% 
    # Preserve sorted order
    mutate(Description = factor(Description, levels = unique(Description)))
  
  # Plot
  df_plot <- ggplot(interesting_ekeggs, aes(x = GeneRatio, y = Description,
                                            color = p.adjust)) +
    geom_point(aes(size = Count)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = title,
      y = "") +
    theme_light() +
    theme(axis.text.y = element_text(size = font_size))
  
  return(df_plot)
  
}



# Searches a given tibble for a term. Basically grep but with fewer characters to type. 

findEGO <- function(ego_tibble, term, print_top_matches = FALSE){
  # Grep search for given term 
  matches <- ego_tibble[grep(term, ego_tibble$Description, ignore.case = TRUE), ]
  
  # Output number of matches
  print(paste0(nrow(matches), " matches found."))
  
  print(matches$Description)
  
  # If TRUE, print top 3 matches along with the genes
  
  if(print_top_matches == TRUE) {
    # If < 4 matches found, print all of them
    if (nrow(matches) < 4) {
      print(paste0("Printing all ", nrow(matches), " matches."))
      
      as.data.frame(matches[ , colnames(matches) %in% c("Description", "GeneRatio", "p.adjust", "geneID")])
    } else {
      
      print(paste0("Matches are: "))
      print(matches$Description)
      
      print(paste0("Printing top 3 matches."))
      
      df_match <- as.data.frame(matches[ , colnames(matches) %in% c("Description", "GeneRatio", "p.adjust", "geneID")])
      
      df_match <- df_match[order(df_match$p.adjust, decreasing = FALSE), ]
      
      head(df_match, n = 3)
    }
  }
  
}




