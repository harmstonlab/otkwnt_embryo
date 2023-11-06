# 03_multiqc

In this folder, we run MultiQC for a quick quality check of all our samples. 

First, we generate an interactive report: 

```bash
multiqc --force -n multiqc_report_interactive \
../00_fastqc/*.zip \
../01_STAR/*Log.final.out \
../02_RSEM/
```

```bash
multiqc --pdf --force --export -n multiqc_report \
../00_fastqc/*.zip \
../01_STAR/*Log.final.out \
../02_RSEM/
```