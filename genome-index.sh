#!/bin/bash

singularity pull docker://evolbioinfo/star:v2.7.6a

#mkdir ref

singularity run star_v2.7.6a.sif --runThreadN 4 \
       	--runMode genomeGenerate \
       	--genomeDir ref/ \
       	--genomeFastaFiles ref.fa



