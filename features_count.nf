#!/usr/bin/env nextflow
process getGenomic_features{

	output:
	file "Homo_sapiens.GRCh38.101.chr.gtf" into FileCount 
	
	script:
	"""
	wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
	gzip -d  Homo_sapiens.GRCh38.101.chr.gtf.gz
	"""
}
 
 process getCount_feature{
	input:
	file "Homo_sapiens.GRCh38.101.chr.gtf" from FileCount
	file "*.bam" from Filebam # à décider
	
	output:
	file "output.counts" into FileCount
	file "output.counts.summary" into FileCount
	
	script
	""" 
 	./subread-v2.0.1.simg -p -t gene -g gene_id -s 0 -a Homo_sapiens.GRCh38.101.chr.gtf -o output.counts *.bam
	"""
}	