#!/bin/bash
#
make -f ./pipeline/makefile \
     UNIPROT=/project/cotton_genomics/uniprot/uniprot_sprot.fasta \
     NAME=B1 \
     MAKER_PROTEIN=/project/cotton_genomics/evidence/Gohir.Gorai.uniport.fa \
     MAKER_ALTEST=/project/cotton_genomics/evidence/gossypium.ests.fa \
     MAKER_MAX_DNA_LEN=1000000 \
     GENOME=/project/cotton_genomics/cotton/genomes/B1/anom_juiced.V1.1.fasta \
     GENEMARK=/project/cotton_genomics/pipeline/gmes_linux_64 \
     MIKADO_SPLIT=30 \
     MIKADO_SCORE=plant.yaml \
     GENOME_SPLIT=13 \
     LIBS="SRR617075 SRR617073 SRR617068 SRR617067 SRR10675236 SRR10675235 SRR10675234 SRR10675237 SRR1174179 SRR617009 SRR617011 SRR617013 SRR959508 SRR2132267 SRR959585 SRR6327757 SRR6327758 SRR6327759 SRR8267554 SRR8267566 SRR8878565 SRR8878526 SRR8878661 SRR8878800 SRR8878534 SRR8878745 SRR8267623 SRR8267616 SRR8267619 SRR8267606 SRR8267582 SRR8267601" \
     SRR617075="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617075.fastq.gz" \
     SRR617073="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617073.fastq.gz" \
     SRR617068="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617068.fastq.gz" \
     SRR617067="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617067.fastq.gz" \
     SRR10675236="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR10675236.fastq.gz" \
     SRR10675235="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR10675235.fastq.gz" \
     SRR10675234="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR10675234.fastq.gz" \
     SRR10675237="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR10675237.fastq.gz" \
     SRR1174179="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR1174179.fastq.gz" \
     SRR617009="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617009.fastq.gz" \
     SRR617011="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617011.fastq.gz" \
     SRR617013="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR617013.fastq.gz" \
     SRR959508="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR959508_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR959508_2.fastq.gz" \
     SRR2132267="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR2132267_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR2132267_2.fastq.gz" \
     SRR959585="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR959585_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR959585_2.fastq.gz" \
     SRR6327757="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR6327757_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR6327757_2.fastq.gz" \
     SRR6327758="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR6327758_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR6327758_2.fastq.gz" \
     SRR6327759="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR6327759_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR6327759_2.fastq.gz" \
     SRR8267554="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267554_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267554_2.fastq.gz" \
     SRR8267566="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267566_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267566_2.fastq.gz" \
     SRR8878565="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878565_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878565_2.fastq.gz" \
     SRR8878526="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878526_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878526_2.fastq.gz" \
     SRR8878661="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878661_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878661_2.fastq.gz" \
     SRR8878800="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878800_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878800_2.fastq.gz" \
     SRR8878534="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878534_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878534_2.fastq.gz" \
     SRR8878745="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878745_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8878745_2.fastq.gz" \
     SRR8267623="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267623_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267623_2.fastq.gz" \
     SRR8267616="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267616_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267616_2.fastq.gz" \
     SRR8267619="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267619_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267619_2.fastq.gz" \
     SRR8267606="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267606_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267606_2.fastq.gz" \
     SRR8267582="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267582_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267582_2.fastq.gz" \
     SRR8267601="/project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267601_1.fastq.gz /project/cotton_genomics/cotton/genomes/B1/RNA-seq/SRR8267601_2.fastq.gz" \
     FILTER=0.37 \
     LINEAGE=eudicots_odb10  queue_extra="-A cotton_genomics" $@
