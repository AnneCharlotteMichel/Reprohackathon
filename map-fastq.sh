#!/bin/bash

singularity pull docker://evolbioinfo/star:v2.7.6a

singularity pull docker://evolbioinfo/samtools:v1.11

singularity run star_v2.7.6a.sif --outSAMstrandField intronMotif \
	--outFilterMismatchNmax 4 \
	--outFilterMultimapNmax 10 \
	--genomeDir ref \
	--readFilesIn <(gunzip -c SRR628582_1.fastq.gz) <(gunzip -c SRR628582_2.fastq.gz) \
	--runThreadN 4 \
	--outSAMunmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outStd BAM_SortedByCoordinate \
	--genomeLoad NoSharedMemory \
	--limitBAMsortRAM 12000000000 \
       	> SRR628582.bam

singularity run samtools_v1.11.sif index *.bam




