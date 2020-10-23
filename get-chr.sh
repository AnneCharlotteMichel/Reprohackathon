#!/bin/bash

singularity pull docker://evolbioinfo/ubuntu:v16.04

liste_chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr19 chr20 chr21 ch22 chrM chrX chrY"

for chromosome in $liste_chromosomes
do 
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/$chromosome.fa.gz
done

gunzip -c *.fa.gz > ref.fa


