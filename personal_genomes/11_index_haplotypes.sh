#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
STAR=/DATA/users/m.magnitov/software/STAR-2.7.10a/bin/Linux_x86_64/STAR

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	# Index genomes for BWA
	bwa index ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.fa
	bwa index ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.fa

	# Index genomes for STAR
	mkdir -p ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_STARgenome
	mkdir -p ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_STARgenome
	${STAR} --runMode genomeGenerate --runThreadN 32 \
		--genomeDir ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_STARgenome \
		--genomeFastaFiles ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.fa \
		--sjdbGTFfile ${DATA_PATH}/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.hap1.gtf
	${STAR} --runMode genomeGenerate --runThreadN 32 \
		--genomeDir ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_STARgenome \
		--genomeFastaFiles ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.fa \
		--sjdbGTFfile ${DATA_PATH}/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.hap2.gtf
done
