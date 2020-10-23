#!/usr/bin/env nextflow
SRAID=SRR628582 


// à tester 
process getFastq{ 
    input : 
    val srr from SRAID

    output :
    file ${srr}_1.fastq.gz into fichiers_fastq
    file ${srr}_2.fastq.gz into fichiers_fastq

    script :
    """
    wget -O ${srr}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${srr}/${srr}.1
    singularity run sratoolkit_v2.5.7.sif fastq-dump --gzip --split-files ./${srr}.sra   
    """

}

process getChr {
    output :
    file ref.fa into genome_file

    script :
    """
    liste_chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr19 chr20 chr21 ch22 chrM chrX chrY"

    for chromosome in $liste_chromosomes
    do 
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/$chromosome.fa.gz
    done

    gunzip -c *.fa.gz > ref.fa
    """
}

process genomeIndex {
    input :
    file ref.fa from genome_file

    output :
    // ? 

    script :
    """
    STAR --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir ref/ \
    --genomeFastaFiles ref.fa
    """
}

process genomeAnnot {
    input:

    output:

    script:
    """
    """
}

process mapFastq {
    input:
    val srr from SRAID

    output:
    file *.bam.bai into mapped_fastq_files

    script:
    """
    STAR --outSAMstrandField intronMotif \ 
        --outFilterMismatchNmax 4 \ 
        --outFilterMultimapNmax 10 \ 
        --genomeDir ref \
        --readFilesIn <(gunzip -c ${srr}_1.fastq.gz) <(gunzip -c ${srr}_2.fastq.gz) \ 
        --runThreadN 4 \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate 
        --outStd BAM_SortedByCoordinate \ 
        --genomeLoad NoSharedMemory \ 
        --limitBAMsortRAM 12000000000 \ 
        > ${srr}.bam
    samtools index *.bam
    """

}