#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
GENCODE_ANNOTATION=/DATA/users/m.magnitov/genomes/gencode.v42.basic.annotation.no_header.gtf

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	paste ${DATA_PATH}/${SAMPLE}/genome/lifted_gencode_annotation_coordinates.hap1.gtf ${GENCODE_ANNOTATION} | grep -v 'error_chrom' |\
		awk '{ print $5"\t"$6"\t"$7"\t"$2"\t"$4 }' > ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap1.left.txt
	paste ${DATA_PATH}/${SAMPLE}/genome/lifted_gencode_annotation_coordinates.hap1.gtf ${GENCODE_ANNOTATION} | grep -v 'error_chrom' |\
		cut -f10- > ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap1.right.txt

	paste ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap1.left.txt \
		${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap1.right.txt >\
		${DATA_PATH}/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.hap1.gtf
        
	paste ${DATA_PATH}/${SAMPLE}/genome/lifted_gencode_annotation_coordinates.hap2.gtf ${GENCODE_ANNOTATION} | grep -v 'error_chrom' |\
		awk '{ print $5"\t"$6"\t"$7"\t"$2"\t"$4 }' > ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap2.left.txt
	paste ${DATA_PATH}/${SAMPLE}/genome/lifted_gencode_annotation_coordinates.hap2.gtf ${GENCODE_ANNOTATION} | grep -v 'error_chrom' |\
		cut -f10- > ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap2.right.txt

	paste ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap2.left.txt \
		${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap2.right.txt >\
		${DATA_PATH}/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.hap2.gtf

	rm ${DATA_PATH}/${SAMPLE}/genome/temp_lifted_gencode_annotation_coordinates.hap*
	rm ${DATA_PATH}/${SAMPLE}/genome/lifted_gencode_annotation_coordinates.hap*
done
