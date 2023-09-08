#!/bin/sh
# Running RSEM

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
  rsem-calculate-expression --paired-end \
                            --alignments 
                            --no-bam-output 
                            -p 20 \
                            /home/qianhui/otkwnt_project/data/data_processed/01_STAR/${i}Aligned.toTranscriptome.out.bam \
                            /home/shared/genomes/dm6/StarIndex/ensembl97/dm6_ensembl97 \
                            $i
done


# Testing code: Run for just 1 condition

#  rsem-calculate-expression --paired-end --alignments --no-bam-output -p 20 /home/qianhui/capstone/X401SC20100446-Z01-F008/analysis/alignReads/O1W4_OE_1Aligned.toTranscriptome.out.bam /home/shared/genomes/dm6/StarIndex/ensembl97/dm6_ensembl97 O1W4_OE_1
