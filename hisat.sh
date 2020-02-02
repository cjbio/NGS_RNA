#!/bin/sh
#RNA-seq workflow with TOPHAT2 and cufflinks.


# - defining folder structure -

#User must specify the path to the RNA-seq project folder
PROJECTFOLDER=/home/chao/rnaseq/ex0072/ex0072-95626531/run1/

#User must specify the name of each sample. i.e. for sample1.fastq.gz and sample2.fastq.gz, write:
SAMPLES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48"

#these folders are auto-generated during running the script. Each stage of the processed data will be saved under one these folders. No user intervention needed.
SCRIPTS=$PROJECTFOLDER/scripts/
FASTQFILES=$PROJECTFOLDER/data/
ALIGNMENTFILES=$PROJECTFOLDER/alignments/
ASSEMBLIESFILES=$PROJECTFOLDER/assemblies/
QUANTIFICATION=$PROJECTFOLDER/quantification
DE=$PROJECTFOLDER/DE/


# - reference genome, index and annotation -

#User to specify the path to these files
GENE_REFERENCE=/lib/GenomeRefs/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/genes.gtf # path to reference annotation gtf or gff file
HiSat_INDEX=$HISAT # path to hisat index
REFERENCEFA=/lib/GenomeRefs/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa # path to reference genome
P=30 # number of threads to use for data analysis (find out number of availible threads on linux using "htop" command)
LIBRARYTYPE=fr-firststrand

# - hisat alignment -

mkdir $SCRIPTS $ALIGNMENTFILES $ASSEMBLIESFILES $QUANTIFICATION $DE

for SAMPLE_ID in $SAMPLES
do
cat > $SCRIPTS/hisat_${SAMPLE_ID}.sh <<EOF
hisat2 --rna-strandness RF --dta -p 30 -x $HISAT -1  <(gzip -cd $FASTQFILES/${SAMPLE_ID}_R1.fastq.gz) -2 <(gzip -cd $FASTQFILES/${SAMPLE_ID}_R2.fastq.gz) -S /dev/stdout 2>> hisatoutput | newsamtools view -bhS -O BAM -@ 2 - > $ALIGNMENTFILES/${SAMPLE_ID}.bam
EOF

bash $SCRIPTS/hisat_${SAMPLE_ID}.sh

done &&

echo "alignment completed"
