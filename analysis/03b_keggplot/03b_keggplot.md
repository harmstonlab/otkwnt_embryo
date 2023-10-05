# kegg_plots
Qian Hui TAN
2023-10-05

- [<span class="toc-section-number">1</span> Pathview
  plots](#pathview-plots)
  - [<span class="toc-section-number">1.1</span> Read in
    data](#read-in-data)
  - [<span class="toc-section-number">1.2</span> Functions](#functions)
- [<span class="toc-section-number">2</span> Prepare
  data](#prepare-data)
- [<span class="toc-section-number">3</span> Plots](#plots)
- [<span class="toc-section-number">4</span> sessionInfo](#sessioninfo)

# Pathview plots

In this notebook, we examine KEGG `pathview` plots for our conditions of
interest.

``` r
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pathview)
})
```

## Read in data

``` r
dds_dm <- readRDS("../output/01_QC/dds_filt.RDS")
ensembl_dm <- readRDS("../output/01_QC/ensembl_genes.RDS")
```

## Functions

``` r
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
  res = results(wald_dds, 
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
    res <- lfcShrink(wald_dds, coef = deseq_coef, 
                                    res = res, 
                                    type = lfcshrinktype, parallel = TRUE)
  }
  
  # Add gene annotations
  res$entrezgene_id = ensembl.genes$entrezgene_id[match(row.names(res), ensembl.genes$gene_id)]
  res$gene_biotype = ensembl.genes$gene_biotype[match(row.names(res), ensembl.genes$gene_id)]
  res$external_gene_name = ensembl.genes$external_gene_name[match(row.names(res), ensembl.genes$gene_id)]
  
  return(res)
}
```

``` r
plot_pathview <- function(res_input, 
                          kegg_pathway, 
                          organism,
                          name){
  
  pathview(gene.data = res_input,
         pathway.id = kegg_pathway,
         kegg.dir = "../figures/03b_keggplot/", 
         species = organism,
         discrete = FALSE, 
         limit = list(gene = 3),
         bins = 20, 
         out.suffix = name, 
         low = "blue", mid = "white", high = "red", 
         na.col = "grey")

}

# Run:
#plot_pathview(rld_c11, kegg_pathway = path_to_plot,
#                organism = "hsa", 
#                name = paste0("c11_", path_name))
```

``` r
# This function takes in the various rlds we've specified, and makes kegg pathview plots for all of them
# note that the rld names and output names have been hard-coded into the function - this isn't meant
# for generic usage, and is customized for this particular dataset because I don't feel like retyping all 
# of this 

plot_multipaths <- function(path_to_plot, path_name) {
  # Plot the various plots
  plot_pathview(res_w4oe_plot, kegg_pathway = path_to_plot,
                organism = "dme", 
                name = paste0("W4OE_", path_name))
  
  plot_pathview(res_o1oe_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O1OE_", path_name))
  
  plot_pathview(res_o1w4oe_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O1W4OE_", path_name))
  
  plot_pathview(res_o2ko_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O2KO_", path_name))
  
  plot_pathview(res_w4ko_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("W4KO_", path_name))
  
  plot_pathview(res_o1ko_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O1KO_", path_name))
  
  plot_pathview(res_o1w4ko_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O1W4KO_", path_name))
  
  plot_pathview(res_o2w4ko_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O2W4KO_", path_name))
  
  plot_pathview(res_o1o2ko_plot, kegg_pathway = path_to_plot,
                organism = "dme", name = paste0("O1O2KO_", path_name))


}
```

# Prepare data

<div class="panel-tabset">

## Prepare for pathview

Here we prepare data for pathview. We color by log2FoldChange because
that makes more sense.

``` r
# Set WT as the comparison
dds <- dds_dm

# Removing lowly expressed genes, only to be done once at the start of the differential expression step
filter = apply(counts(dds, normalized = TRUE), 1, function(x){ mean(x) >= 10 })
dds = dds[filter, ]

# Rerun DESeq
dds <- DESeq(dds, test = "Wald", parallel = TRUE)
```

    using pre-existing size factors

    estimating dispersions

    gene-wise dispersion estimates: 6 workers

    mean-dispersion relationship

    final dispersion estimates, fitting model and testing: 6 workers

``` r
# Sanity check - all rownames and sampleid in colData must match. 
# if this is FALSE, check the order of data_mat and experimental_metadata - 
# something may have gone wrong. 
all(rownames(colData(dds)) == colData(dds)$sample_id)
```

    [1] TRUE

``` r
resultsNames(dds)
```

     [1] "Intercept"                 "condition_O1_KO_vs_Ctrl"  
     [3] "condition_O1_OE_vs_Ctrl"   "condition_O1O2_KO_vs_Ctrl"
     [5] "condition_O1W4_KO_vs_Ctrl" "condition_O1W4_OE_vs_Ctrl"
     [7] "condition_O2_KO_vs_Ctrl"   "condition_O2W4_KO_vs_Ctrl"
     [9] "condition_W4_KO_vs_Ctrl"   "condition_W4_OE_vs_Ctrl"  

A sanity check: rownames and sampleid must match

## W4OE

``` r
res_w4oe <- get_dds_res(dds,
                      contrast = c("condition", "W4_OE", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_W4_OE_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_w4oe)
```

    log2 fold change (MAP): condition W4_OE vs Ctrl 
    Wald test p-value: condition W4 OE vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE    pvalue      padj
                <numeric>      <numeric> <numeric> <numeric> <numeric>
    FBgn0000003   277.736   -0.000588810 0.0238343  0.643645  0.999566
    FBgn0000008  2202.628   -0.003487071 0.0238111  0.427183  0.999566
    FBgn0000014  4638.107   -0.002848772 0.0238974  0.359940  0.999566
    FBgn0000015  2851.378   -0.001180215 0.0236399  0.262677  0.999566
    FBgn0000017  8384.411   -0.004489117 0.0241394  0.254812  0.999566
    FBgn0000018   264.936   -0.000987052 0.0236982  0.739866  0.999566
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
summary(res_w4oe$log2FoldChange)
```

         Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    -9.076734 -0.001473 -0.000146  0.001389  0.001370  9.801361 

Convert rownames to entrez:

``` r
res_w4oe_plot <- as.matrix(res_w4oe["log2FoldChange"])
rownames(res_w4oe_plot) <- res_w4oe$entrezgene_id

head(res_w4oe_plot)
```

            log2FoldChange
    3771948  -0.0005888097
    43852    -0.0034870712
    42037    -0.0028487717
    47763    -0.0011802153
    45821    -0.0044891172
    44793    -0.0009870516

## O1OE

``` r
res_o1oe <- get_dds_res(dds,
                      contrast = c("condition", "O1_OE", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O1_OE_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o1oe)
```

    log2 fold change (MAP): condition O1_OE vs Ctrl 
    Wald test p-value: condition O1 OE vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE    pvalue      padj
                <numeric>      <numeric> <numeric> <numeric> <numeric>
    FBgn0000003   277.736     -0.0452588  0.201101  0.590052  0.724084
    FBgn0000008  2202.628     -0.1128156  0.114769  0.249971  0.428701
    FBgn0000014  4638.107     -0.1437277  0.154156  0.210739  0.390395
    FBgn0000015  2851.378     -0.1779211  0.146659  0.127490  0.295046
    FBgn0000017  8384.411      0.1362630  0.127616  0.194355  0.371704
    FBgn0000018   264.936      0.1260715  0.155487  0.267586  0.446323
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

Convert rownames to entrez:

``` r
res_o1oe_plot <- as.matrix(res_o1oe["log2FoldChange"])
rownames(res_o1oe_plot) <- res_o1oe$entrezgene_id

head(res_o1oe_plot)
```

            log2FoldChange
    3771948    -0.04525882
    43852      -0.11281560
    42037      -0.14372774
    47763      -0.17792113
    45821       0.13626304
    44793       0.12607150

## O1W4_OE

``` r
res_o1w4oe <- get_dds_res(dds,
                      contrast = c("condition", "O1W4_OE", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O1W4_OE_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o1w4oe)
```

    log2 fold change (MAP): condition O1W4_OE vs Ctrl 
    Wald test p-value: condition O1W4 OE vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE    pvalue      padj
                <numeric>      <numeric> <numeric> <numeric> <numeric>
    FBgn0000003   277.736      0.0439553  0.189300 0.5607929  0.773137
    FBgn0000008  2202.628     -0.1632005  0.118686 0.0943964  0.301360
    FBgn0000014  4638.107     -0.2006758  0.165940 0.0805362  0.278760
    FBgn0000015  2851.378     -0.2398033  0.156198 0.0450289  0.208780
    FBgn0000017  8384.411     -0.0798753  0.120242 0.4126817  0.663494
    FBgn0000018   264.936      0.0639870  0.141612 0.5322306  0.754429
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_o1w4oe_plot <- as.matrix(res_o1w4oe["log2FoldChange"])
rownames(res_o1w4oe_plot) <- res_o1w4oe$entrezgene_id

head(res_o1w4oe_plot)
```

            log2FoldChange
    3771948     0.04395535
    43852      -0.16320053
    42037      -0.20067578
    47763      -0.23980335
    45821      -0.07987525
    44793       0.06398699

## O2_KO

``` r
res_o2ko <- get_dds_res(dds,
                      contrast = c("condition", "O2_KO", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O2_KO_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o2ko)
```

    log2 fold change (MAP): condition O2_KO vs Ctrl 
    Wald test p-value: condition O2 KO vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE    pvalue      padj
                <numeric>      <numeric> <numeric> <numeric> <numeric>
    FBgn0000003   277.736      0.0354818  0.163585 0.5346166  0.743864
    FBgn0000008  2202.628     -0.2809791  0.130232 0.0061172  0.104725
    FBgn0000014  4638.107     -0.1680978  0.160736 0.0990122  0.343874
    FBgn0000015  2851.378     -0.2962987  0.167655 0.0124258  0.144557
    FBgn0000017  8384.411     -0.0339997  0.110225 0.6831389  0.838432
    FBgn0000018   264.936     -0.0283121  0.127792 0.7366822  0.870658
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_o2ko_plot <- as.matrix(res_o2ko["log2FoldChange"])
rownames(res_o2ko_plot) <- res_o2ko$entrezgene_id

head(res_o2ko_plot)
```

            log2FoldChange
    3771948     0.03548184
    43852      -0.28097910
    42037      -0.16809782
    47763      -0.29629871
    45821      -0.03399967
    44793      -0.02831211

## W4KO

``` r
res_w4ko <- get_dds_res(dds,
                      contrast = c("condition", "W4_KO", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_W4_KO_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_w4ko)
```

    log2 fold change (MAP): condition W4_KO vs Ctrl 
    Wald test p-value: condition W4 KO vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE      pvalue        padj
                <numeric>      <numeric> <numeric>   <numeric>   <numeric>
    FBgn0000003   277.736      -0.141803  0.354595 5.90531e-01 6.70096e-01
    FBgn0000008  2202.628      -0.497346  0.127125 5.09849e-05 2.52010e-04
    FBgn0000014  4638.107      -0.673087  0.179521 6.06008e-05 2.90493e-04
    FBgn0000015  2851.378      -0.558772  0.160144 2.49977e-04 9.61019e-04
    FBgn0000017  8384.411       0.749580  0.142089 4.16941e-08 7.22676e-07
    FBgn0000018   264.936       0.330460  0.181212 5.06965e-02 8.77894e-02
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_w4ko_plot <- as.matrix(res_w4ko["log2FoldChange"])
rownames(res_w4ko_plot) <- res_w4ko$entrezgene_id

head(res_w4ko_plot)
```

            log2FoldChange
    3771948     -0.1418025
    43852       -0.4973455
    42037       -0.6730870
    47763       -0.5587718
    45821        0.7495798
    44793        0.3304602

## O1KO

``` r
res_o1ko <- get_dds_res(dds,
                      contrast = c("condition", "O1_KO", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O1_KO_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o1ko)
```

    log2 fold change (MAP): condition O1_KO vs Ctrl 
    Wald test p-value: condition O1 KO vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE      pvalue       padj
                <numeric>      <numeric> <numeric>   <numeric>  <numeric>
    FBgn0000003   277.736      -0.147006  0.324819 0.486242217 0.60396772
    FBgn0000008  2202.628      -0.438989  0.126846 0.000288402 0.00244429
    FBgn0000014  4638.107      -0.493094  0.178416 0.002329831 0.01117296
    FBgn0000015  2851.378      -0.307247  0.155907 0.035258763 0.08100575
    FBgn0000017  8384.411       0.136917  0.135543 0.279865900 0.39835319
    FBgn0000018   264.936       0.275995  0.178041 0.084834692 0.15784351
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_o1ko_plot <- as.matrix(res_o1ko["log2FoldChange"])
rownames(res_o1ko_plot) <- res_o1ko$entrezgene_id

head(res_o1ko_plot)
```

            log2FoldChange
    3771948     -0.1470061
    43852       -0.4389889
    42037       -0.4930944
    47763       -0.3072466
    45821        0.1369174
    44793        0.2759952

## O1W4KO

``` r
res_o1w4ko <- get_dds_res(dds,
                      contrast = c("condition", "O1W4_KO", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O1W4_KO_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o1w4ko)
```

    log2 fold change (MAP): condition O1W4_KO vs Ctrl 
    Wald test p-value: condition O1W4 KO vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE      pvalue        padj
                <numeric>      <numeric> <numeric>   <numeric>   <numeric>
    FBgn0000003   277.736       0.325157  0.389435 2.69614e-01 3.47461e-01
    FBgn0000008  2202.628      -0.866411  0.129076 5.52850e-12 1.26224e-10
    FBgn0000014  4638.107      -1.059318  0.181496 7.26957e-10 9.36582e-09
    FBgn0000015  2851.378      -0.847002  0.162031 6.17003e-08 4.61510e-07
    FBgn0000017  8384.411       0.729191  0.141849 1.01912e-07 7.22523e-07
    FBgn0000018   264.936       0.214726  0.181268 2.09404e-01 2.80439e-01
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_o1w4ko_plot <- as.matrix(res_o1w4ko["log2FoldChange"])
rownames(res_o1w4ko_plot) <- res_o1w4ko$entrezgene_id

head(res_o1w4ko_plot)
```

            log2FoldChange
    3771948      0.3251566
    43852       -0.8664106
    42037       -1.0593184
    47763       -0.8470019
    45821        0.7291907
    44793        0.2147259

## O2W4KO

``` r
res_o2w4ko <- get_dds_res(dds,
                      contrast = c("condition", "O2W4_KO", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O2W4_KO_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o2w4ko)
```

    log2 fold change (MAP): condition O2W4_KO vs Ctrl 
    Wald test p-value: condition O2W4 KO vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE      pvalue        padj
                <numeric>      <numeric> <numeric>   <numeric>   <numeric>
    FBgn0000003   277.736      0.0131126  0.334764 9.60842e-01 9.74027e-01
    FBgn0000008  2202.628     -0.4224889  0.126226 4.32513e-04 1.85113e-03
    FBgn0000014  4638.107     -0.4945494  0.177252 2.55924e-03 7.89228e-03
    FBgn0000015  2851.378     -0.4347631  0.158417 3.74750e-03 1.08202e-02
    FBgn0000017  8384.411      0.6870706  0.141992 4.23215e-07 7.86379e-06
    FBgn0000018   264.936      0.1308322  0.175861 4.22738e-01 5.22620e-01
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_o2w4ko_plot <- as.matrix(res_o2w4ko["log2FoldChange"])
rownames(res_o2w4ko_plot) <- res_o2w4ko$entrezgene_id

head(res_o2w4ko_plot)
```

            log2FoldChange
    3771948     0.01311263
    43852      -0.42248893
    42037      -0.49454939
    47763      -0.43476314
    45821       0.68707058
    44793       0.13083225

## O1O2KO

``` r
res_o1o2ko <- get_dds_res(dds,
                      contrast = c("condition", "O1O2_KO", "Ctrl"), 
                      ensembl.genes = ensembl_dm)
```

    [1] "condition_O1O2_KO_vs_Ctrl"

    using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        sequence count data: removing the noise and preserving large differences.
        Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
head(res_o1o2ko)
```

    log2 fold change (MAP): condition O1O2_KO vs Ctrl 
    Wald test p-value: condition O1O2 KO vs Ctrl 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE      pvalue       padj
                <numeric>      <numeric> <numeric>   <numeric>  <numeric>
    FBgn0000003   277.736     0.09509622  0.335533 0.708975049 0.78790514
    FBgn0000008  2202.628    -0.31906299  0.125323 0.007244987 0.02292600
    FBgn0000014  4638.107    -0.62445895  0.179465 0.000167039 0.00122963
    FBgn0000015  2851.378    -0.44272600  0.158685 0.003091834 0.01186178
    FBgn0000017  8384.411     0.19648086  0.137415 0.131261535 0.21858641
    FBgn0000018   264.936     0.00250876  0.175559 0.985355926 0.99060002
                entrezgene_id   gene_biotype external_gene_name
                    <integer>    <character>        <character>
    FBgn0000003       3771948          ncRNA     7SLRNA:CR32864
    FBgn0000008         43852 protein_coding                  a
    FBgn0000014         42037 protein_coding              abd-A
    FBgn0000015         47763 protein_coding              Abd-B
    FBgn0000017         45821 protein_coding                Abl
    FBgn0000018         44793 protein_coding                abo

``` r
res_o1o2ko_plot <- as.matrix(res_o1o2ko["log2FoldChange"])
rownames(res_o1o2ko_plot) <- res_o1o2ko$entrezgene_id

head(res_o1o2ko_plot)
```

            log2FoldChange
    3771948     0.09509622
    43852      -0.31906299
    42037      -0.62445895
    47763      -0.44272600
    45821       0.19648086
    44793       0.00250876

</div>

# Plots

I’d like the figures to be in a separate kegg directory, but somehow
that option (`kegg.dir`) is not working and insists on saving to the
default working directory. For now we just plot everything and then use
the terminal to shift all \*.png files into a separate figures directory
later.

Not an ideal solution, but this will do for now.

``` r
plot_multipaths("04391", "hippo")
```

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.W4OE_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O1OE_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O1W4OE_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O2KO_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.W4KO_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O1KO_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O1W4KO_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O2W4KO_hippo.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04391.O1O2KO_hippo.png

``` r
plot_multipaths("04013", "MAPK")
```

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.W4OE_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O1OE_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O1W4OE_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O2KO_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.W4KO_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O1KO_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O1W4KO_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O2W4KO_MAPK.png

    Info: some node width is different from others, and hence adjusted!

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04013.O1O2KO_MAPK.png

    Info: some node width is different from others, and hence adjusted!

``` r
plot_multipaths("04330", "Notch")
```

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.W4OE_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O1OE_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O1W4OE_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O2KO_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.W4KO_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O1KO_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O1W4KO_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O2W4KO_Notch.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04330.O1O2KO_Notch.png

``` r
plot_multipaths("04350", "TGFbeta")
```

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.W4OE_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O1OE_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O1W4OE_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O2KO_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.W4KO_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O1KO_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O1W4KO_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O2W4KO_TGFbeta.png

    Info: Getting gene ID data from KEGG...

    Info: Done with data retrieval!

    Info: Working in directory /Users/qh_tan/qianhui/work/otkwnt_embryo/analysis/03b_keggplot

    Info: Writing image file dme04350.O1O2KO_TGFbeta.png

# sessionInfo

``` r
sessionInfo()
```

    R version 4.2.2 (2022-10-31)
    Platform: aarch64-apple-darwin20 (64-bit)
    Running under: macOS Ventura 13.1

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
     [1] pathview_1.38.0             DESeq2_1.38.3              
     [3] SummarizedExperiment_1.28.0 Biobase_2.58.0             
     [5] MatrixGenerics_1.10.0       matrixStats_1.0.0          
     [7] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
     [9] IRanges_2.32.0              S4Vectors_0.36.2           
    [11] BiocGenerics_0.44.0         forcats_1.0.0              
    [13] stringr_1.5.0               dplyr_1.1.3                
    [15] purrr_1.0.1                 readr_2.1.3                
    [17] tidyr_1.3.0                 tibble_3.2.1               
    [19] ggplot2_3.4.2               tidyverse_1.3.2            

    loaded via a namespace (and not attached):
     [1] googledrive_2.0.0      colorspace_2.1-0       XVector_0.38.0        
     [4] fs_1.6.2               rstudioapi_0.14        bit64_4.0.5           
     [7] AnnotationDbi_1.60.2   fansi_1.0.4            mvtnorm_1.1-3         
    [10] apeglm_1.20.0          lubridate_1.9.1        xml2_1.3.5            
    [13] codetools_0.2-19       cachem_1.0.8           geneplotter_1.76.0    
    [16] knitr_1.42             jsonlite_1.8.7         broom_1.0.3           
    [19] annotate_1.76.0        dbplyr_2.3.3           png_0.1-8             
    [22] graph_1.76.0           compiler_4.2.2         httr_1.4.6            
    [25] backports_1.4.1        Matrix_1.5-4.1         fastmap_1.1.1         
    [28] gargle_1.3.0           cli_3.6.1              htmltools_0.5.4       
    [31] tools_4.2.2            coda_0.19-4            gtable_0.3.3          
    [34] glue_1.6.2             GenomeInfoDbData_1.2.9 Rcpp_1.0.11           
    [37] bbmle_1.0.25           cellranger_1.1.0       vctrs_0.6.3           
    [40] Biostrings_2.66.0      xfun_0.37              rvest_1.0.3           
    [43] timechange_0.2.0       lifecycle_1.0.3        XML_3.99-0.14         
    [46] googlesheets4_1.0.1    org.Hs.eg.db_3.16.0    zlibbioc_1.44.0       
    [49] MASS_7.3-58.2          scales_1.2.1           hms_1.1.3             
    [52] parallel_4.2.2         KEGGgraph_1.58.3       RColorBrewer_1.1-3    
    [55] curl_5.0.1             yaml_2.3.7             memoise_2.0.1         
    [58] emdbook_1.3.12         bdsmatrix_1.3-6        stringi_1.7.12        
    [61] RSQLite_2.3.1          BiocParallel_1.32.6    rlang_1.1.1           
    [64] pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.20         
    [67] lattice_0.20-45        bit_4.0.5              tidyselect_1.2.0      
    [70] plyr_1.8.8             magrittr_2.0.3         R6_2.5.1              
    [73] generics_0.1.3         DelayedArray_0.24.0    DBI_1.1.3             
    [76] pillar_1.9.0           haven_2.5.1            withr_2.5.0           
    [79] KEGGREST_1.38.0        RCurl_1.98-1.12        modelr_0.1.10         
    [82] crayon_1.5.2           utf8_1.2.3             tzdb_0.3.0            
    [85] rmarkdown_2.20         locfit_1.5-9.7         grid_4.2.2            
    [88] readxl_1.4.1           blob_1.2.4             Rgraphviz_2.42.0      
    [91] reprex_2.0.2           digest_0.6.33          xtable_1.8-4          
    [94] numDeriv_2016.8-1.1    munsell_0.5.0         
