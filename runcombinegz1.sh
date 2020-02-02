#!/bin/sh

SAMPLES="s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s30 s31 s32 s33 s34 s35 s36 s37 s38 s39 s40 s41 s42 s43 s44 s45 s46 s47"

for SAMPLES_ID in $SAMPLES
do

cat ${SAMPLES_ID}_R1.fastq.gz ${SAMPLES_ID}_R1_2.fastq.gz > all_${SAMPLES_ID}_R1.fastq.gz
cat ${SAMPLES_ID}_R2.fastq.gz ${SAMPLES_ID}_R2_2.fastq.gz > all_${SAMPLES_ID}_R2.fastq.gz

#cat DMEM-${SAMPLES_ID}-p16_S*_L001_R1_001.fastq.gz DMEM-${SAMPLES_ID}-p16_S*_L002_R1_001.fastq.gz DMEM-${SAMPLES_ID}-p16_S*_L003_R1_001.fastq.gz DMEM-${SAMPLES_ID}-p16_S*_L004_R1_001.fastq.gz > DMEM_${SAMPLES_ID}_p16_R1.fastq.gz
#cat DMEM-${SAMPLES_ID}-p16_S*_L001_R2_001.fastq.gz DMEM-${SAMPLES_ID}-p16_S*_L002_R2_001.fastq.gz DMEM-${SAMPLES_ID}-p16_S*_L003_R2_001.fastq.gz DMEM-${SAMPLES_ID}-p16_S*_L004_R2_001.fastq.gz > DMEM_${SAMPLES_ID}_p16_R2.fastq.gz 


done

