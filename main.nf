#!/usr/bin/env nextflow

// Liste des identifiants (SRAID) des fichiers à télécharger
SRAID = Channel.of("SRR628582","SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")

// Liste des chromosomes à télécharger
liste_chromosomes = Channel.of("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")

// Liste des fichiers bam pour le process deseq
SRAID2 = Channel.of('c("SRR628582.bam","SRR628584.bam", "SRR628586.bam", "SRR628588.bam", "SRR628589.bam", "SRR628583.bam", "SRR628585.bam", "SRR628587.bam")')

// Mise en paramètre du fichier coldata.csv (fichier renseignant sur le plan d'expérience) lors du run 
params.design_file="/mnt/Reprohackathon/coldata.csv"
fichier_mutant = Channel.of(params.design_file)


// Téléchargement des fichiers sra à partir de la liste des SRAID
process getFastq{
    input :
    val srr from SRAID
 
    output :
    tuple val(srr), file("*.sra") into sra_files
 
    script :
    """
    wget -O ${srr}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${srr}/${srr}.1
    """
}


// Conversion des fichiers SRA en fichier FATSQ compressé
process dumpFastq{
    publishDir "fastq/"
    input :
    tuple val(srr), file(sra) from sra_files
 
    output :
    tuple val(srr), file("*_1.fastq.gz"), file("*_2.fastq.gz") into fastq_files
 
    script :
    """
    fastq-dump --gzip --split-files ${srr}.sra  
    """
}


// Téléchargement et décompression des séquences des chromosomes humains à partir de la liste des chromosomes
process getChr {
    input :
    val chr from liste_chromosomes
 
    output :
    file "*.fa" into genome_files
 
    script :
    """
    wget -O chr${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    gunzip chr${chr}.fa.gz
    """
}


// Concaténation des séquences des chromosomes en un seul fichier
process concatChr {
    input:
    file gen from genome_files.collect()
    
    output:
    file "ref.fa" into genome_concat
    
    script :
    """
    cat *.fa > ref.fa
    """
}
 

// Indexation du génome humain à partir du fichier avec les séquences des chromosomes concaténées
process genomeIndex {
    input :
    file gen from genome_concat
 
    output :
    file "ref" into star_index
 
    script :
    """
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir ref \
    --genomeFastaFiles ${gen}
    """
}
 

// Mapping des fichiers FASTQ compressés avec le génome index (alignement des séquences avec le génôme humain) 
process mapFastq {
    publishDir "bam/"

    input:
    tuple val(srr), file(f1), file(f2) from fastq_files
    path ref from star_index
 
    output:
    file "${srr}.bam" into mapped_fastq_files_1, mapped_fastq_files_2
 
    script:
    """
    STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ${ref} \
        --readFilesIn <(gunzip -c ${f1}) <(gunzip -c ${f2}) \
        --runThreadN ${task.cpus} \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM ${task.memory.toBytes()} \
	 >${srr}.bam    
    """
}


// Indexation des fichiers bam
process samFastq {
    publishDir "bai/"

    input:
    file bam from mapped_fastq_files_1
 
    output:
    file "*.bai" into sam_fastq_files
 
    script:
    """
    samtools index *.bam
    """
}
 

// Téléchargement des annotations du génome
process getGenomic_features{
    output:
    file "*.gtf" into annotation 
    
    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gzip -d  Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}


// Comptage des reads
process getCount_feature{
    publishDir "feat_output/"
    input:
    file gtf from annotation
    file bam from mapped_fastq_files_2.collect()
    file bai from sam_fastq_files.collect()
    
    output:
    file "output.counts" into FileCount
    file "output.counts.summary" into logsFileCount
    
    script:
    """ 
     featureCounts -p -t gene -g gene_id -s 0 -a ${gtf} -o output.counts ${bam}
    """
}  


// Analyse statistique à l'aide de R et de la librarie DESeq2
process deseq {
        publishDir "deseq_resultat/"

        input :
        val srr from SRAID2
        path counts from FileCount
        path data from fichier_mutant

        output :
        file "deseq_result.csv" into result_deseq

        script :
        """
        #!/usr/bin/env Rscript
        library("DESeq2")
        cts = read.table("${counts}", sep="\t", header=TRUE, row.names="Geneid")
        cts <- cts[,${srr}]
        coldata = read.csv("${data}", sep=",",header=TRUE, row.names="seq_name")
        dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]
        dds <- DESeq(dds)
        res <- results(dds)
	jpeg(file="saving_plot1.jpeg")
        plotMA(res, ylim=c(-2,2))
	dev.off()
	resOrdered <- res[order(res$padj),]
	write.csv(as.data.frame(resOrdered), file="deseq_result.csv")
	res_filt <- subset(resOrdered, padj < 0.01)
	write.csv(as.data.frame(res_filt), file="deseq_filt.csv")
         """
}


// Sélection des gènes ayant une p-valeur ajustée inférieure à 0.01 (et donc présentant une différence d'expression significative entre les cancers de classe 1 et 2) dans un fichier csv
process filter{
        publishDir "filter_result/"

        input:
        file result from result_deseq
        
        output :
        file "deseq_filtered_res.csv"

        script:
        """
        #!/usr/bin/env Rscript
        library("tidyverse")
        res = read.csv("$result", sep=",", header=TRUE)
        res = filter(res,padj<0.01)
        write.csv(res, file="deseq_filtered_res.csv")
        """
}
