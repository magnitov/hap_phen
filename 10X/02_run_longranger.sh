#!/bin/bash
INPUT_PATH=/DATA/users/m.magnitov/hap_phen/10X
VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0
LONGRANGER=/DATA/users/m.magnitov/software/longranger-2.2.2/longranger
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar
SAMTOOLS=/DATA/users/m.magnitov/software/samtools-1.15/samtools

###
### Run longranger pipeline
###
# LCLs from the lab
for SAMPLE in NA18983 HG01241 HG02601 HG03464
do
	${LONGRANGER} wgs --sample ${SAMPLE} --id ${SAMPLE} --fastqs ${INPUT_PATH}/fastq/${SAMPLE} --reference ${GENOME} --precalled ${VCF_PATH}/${SAMPLE}/${SAMPLE}.vcf --localcores 32 --localmem 256
done

# LCLs from public data
for SAMPLE in NA12878 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	${LONGRANGER} wgs --sample bamtofastq --id ${SAMPLE} --fastqs ${INPUT_PATH}/fastq/${SAMPLE} --reference ${GENOME} --precalled ${VCF_PATH}/${SAMPLE}/${SAMPLE}.vcf --localcores 32 --localmem 256
done

###
### Remove duplicates and reads without barcode from longranger BAM files
###
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	java -jar ${PICARD} MarkDuplicates \
		INPUT=${INPUT_PATH}/${SAMPLE}/outs/phased_possorted_bam.bam \
		OUTPUT=${INPUT_PATH}/${SAMPLE}/outs/phased_possorted_bam.rmdup.bam \
		METRICS_FILE=${INPUT_PATH}/${SAMPLE}/outs/picard.dedup.metrics.txt \
		REMOVE_DUPLICATES=true ASSUME_SORTED=true

	${SAMTOOLS} view -@ 64 -h -S -b -d BX -o ${INPUT_PATH}/${SAMPLE}/outs/phased_possorted_bam.rmdup.rmbx.bam ${INPUT_PATH}/${SAMPLE}/outs/phased_possorted_bam.rmdup.bam     
	${SAMTOOLS} index -@ 64 ${INPUT_PATH}/${SAMPLE}/outs/phased_possorted_bam.rmdup.rmbx.bam
	rm ${INPUT_PATH}/${SAMPLE}/outs/phased_possorted_bam.rmdup.bam
done

