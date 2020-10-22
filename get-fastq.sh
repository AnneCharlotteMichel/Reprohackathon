#!/bin/bash

singularity pull docker://evolbioinfo/sratoolkit:v2.5.7

SRAID=SRR628582 

wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1


singularity run sratoolkit_v2.5.7.sif fastq-dump --gzip --split-files ./SRR628582.sra                                               #permet de télécharger le fichier fastq SRR628582 (utilisé pour le test) à partir du Singularity.
