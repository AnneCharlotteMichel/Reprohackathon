#!/bin/bash

wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
gzip -d  Homo_sapiens.GRCh38.101.chr.gtf.gz

singularity pull docker://evolbioinfo/subread:v2.0.1  
singularity run subread_v2.0.1.sif -p -t gene -g gene_id -s 0 -a Homo_sapiens.GRCh38.101.chr.gtf -o output.counts *.bam
