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
    liste_chromosomes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 20 21 22 MT X Y"

    for chromosome in $liste_chromosomes
    do
            wget chr$chromosome.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$chromosome.fa.gz
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