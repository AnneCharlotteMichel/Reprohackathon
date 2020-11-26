# Reprohackathon
## Présentation du projet
Ce repertoire Github contient le travail du Groupe 7 pour le projet Reprohackathon dans le cadre du M2 AMI2B. Le but du projet est de reproduire certaines parties d'un cas d'analyse de RNAseq (Harbour et al. (Nat. Genet. 2013), Furney et al. (Cancer Discov. 2013)) et d'utiliser un système de management des workflow (Nextflow), des containers (Singularity) et Git. 

Pour cela nous allons analyser le jeu de données utilisé dans ces deux articles et essayer de trouver des gènes dont l'expression est différente. 


## Description du workflow
Notre workflow se trouve dans le fichier main.nf, il s'accompagne du fichier nextflow.config et du fichier coldata.csv. 
* Le process *getFastq* permet de télécharger les fichiers SRA depuis le NCBI. 
* Le process *dumpFastq* permet de convertir les fichiers SRA en FASTQ.
* Le process *getChr* permet de télécharger depuis Ensembl les séquences des chromosomes humains nécessaire pour générer un index sur le génome.
* Le process *concatChr* permet de décompresser les séquences des chromosomes et de les mettre dans un unique fichier.
* Le process *genomeIndex* crée un index sur le génome à l'aide de STAR.
* Le process *mapFastq* indexe les fichiers FATSTQ sur le génome.
* Le process *samFastq* permet de changer le format des alignement avec SAMTOOLS. 
* Le process *getGenomic_features* permet de télécharger les fichiers d'annotations du génome. 
* Le process *getCount_feature* permet de compter les reads alignés sur chaque gène.
* Le process *deseq* permet de faire des tests statistiques pour savoir si une séquences est plus représentée qu'une autre.

Le fichier coldata.csv permet de savoir si le gène SF2B1 est muté selon le fichier SRA considéré. 


## Utilisation sur une machine virtuelle (VM)
* Lancer une VM avec 64Go de RAM et 16 coeurs depuis Biosphère (IFB)
* Activer conda avec la commande suivante :
```
conda activate
```
* Installer Singularity s'il n'est pas déjà installé :
```
conda install singularity=3.6.3
```
* Installer R
```
sudo apt install r-base
```
* Se placer dans au niveau de /mnt/
```
cd /mnt/
```
* Cloner ce répertoire:
```
git clone https://github.com/AnneCharlotteMichel/Reprohackathon
```
* Se placer dans le répertoire Reprohackathon
```
cd Reprohackathon
```
* Autoriser l'exécution des scripts main.nf et nextflow.config
```
chmod +x main.nf
chmod +x nextflow.config
```
* Lancer le script main.nf
```
nextflow run main.nf
```
