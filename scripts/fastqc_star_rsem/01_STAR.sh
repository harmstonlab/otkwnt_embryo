#!/bin/sh
# Running STAR

for i in Ctrl_1 Ctrl_2 Ctrl_3 \
O1_KO_1 O1_KO_2 O1_KO_3 \
O1_OE_1 O1_OE_2 O1_OE_3 \
W4_KO_1 W4_KO_2 W4_KO_3 \
W4_OE_1 W4_OE_2 W4_OE_3 \
O1W4_KO_1 O1W4_KO_2 O1W4_KO_3\
O1W4_OE_1 O1W4_OE_2 O1W4_OE_3 \
O2W4_KO_1 O2W4_KO_2 O2W4_KO_3 \
O2_KO_1 O2_KO_2 O2_KO_3 \
O1O2_KO_1 O1O2_KO_2 O1O2_KO_3 \

do
  STAR  --genomeDir /home/shared/genomes/dm6/StarIndex/ensembl97/ \
        --sjdbGTFfile /home/shared/genomes/dm6/StarIndex/ensembl97/Drosophila_melanogaster.BDGP6.22.97.chr.gtf \
        --readFilesIn /home/qianhui/otkwnt_project/data/data_raw/X401SC20100446-Z01-F008/raw_data/${i}_1.fq.gz \
        /home/qianhui/otkwnt_project/data/data_raw/X401SC20100446-Z01-F008/raw_data/${i}_2.fq.gz  \
        --runThreadN 20 \
        --twopassMode Basic \
        --outWigType bedGraph \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --readFilesCommand zcat \
        --runDirPerm All_RWX \
        --outFileNamePrefix /home/qianhui/otkwnt_project/data/data_processed/01_STAR/${i}
done


# Testing code: Run for just 1 condition and check output

# STAR --genomeDir /home/shared/genomes/dm6/StarIndex/ensembl97/ --sjdbGTFfile /home/shared/genomes/dm6/StarIndex/ensembl97/Drosophila_melanogaster.BDGP6.22.97.chr.gtf --readFilesIn O1W4_OE_1_1.fq.gz O1W4_OE_1_2.fq.gz  --runThreadN 10 --twopassMode Basic --outWigType bedGraph --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --readFilesCommand zcat --runDirPerm All_RWX --outFileNamePrefix O1W4_OE_1
