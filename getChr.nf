#!/usr/bin/env nextflow
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