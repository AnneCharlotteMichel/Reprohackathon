#!/usr/bin/env nextflow

SRAID = Channel.of("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
 
liste_chromosomes = Channel.of("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")
 
process getFastq{
    input :
    val srr from SRAID
 
    output :
    tuple var(srr), file("*_1.fastq.gz"), file("*_2.fastq.gz") into fastq_files
 
    script :
    """
    wget -O ${srr}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${srr}/${srr}.1
    fastq-dump --gzip --split-files ./${srr}.sra  
    """
 
}
 
process getChr {
    input :
    val chr from liste_chromosome
 
    output :
    file "*.fa.gz" into genome_files
 
    script :
    """
    wget chr$chromosome.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$chromosome.fa.gz
    """
}
 
process concatChr{
    input :
    val chr from genome_files.toList()
 
    output:
    file "ref.fa" into genome

    script :
    """
    gunzip -c *.fa.gz > ref.fa
    """
}
 
process genomeIndex {
    input :
    file gen from genome
 
    output :
    file "ref" into star_index
 
    script :
    """
    STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ref \
    --genomeFastaFiles ${gen}
    """
}
 
 
process mapFastq {
    input:
    tuple var(srr), file(f1), file(f2) from fastq_files
 
    output:
    file bam into mapped_fastq_files
 
    script:
    """
    STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ref \
        --readFilesIn <(gunzip -c ${f1}) <(gunzip -c ${f2}) \
        --runThreadN 2 \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 12000000000 \
        > ${srr}.bam
    samtools index *.bam
    """
}
 
getUrl = Channel.of("ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz")
 
process getGenomic_features{
    input:
    val url from getUrl
 
    output:
    file gtf into annotation 
    
    script:
    """
    wget     gzip -d  Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}
 
 process getCount_feature{
    input:
    file gtf from annotation
    file bam from Filebam.toList()
    
    output:
    file "output.counts" into FileCount
    file "output.counts.summary" into logsFileCount
    
    script
    """ 
     featureCounts -p -t gene -g gene_id -s 0 -a gtf -o output.counts bam
    """
}    