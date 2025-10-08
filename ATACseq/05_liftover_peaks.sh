# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
GENOMES_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes

# Merge all peaks from five cell lines into unified peak set
cat ${DATA_PATH}/peaks/*_peaks.canonical.replicated.no_blacklist.bed | sort -k1,1 -k2,2n > ${DATA_PATH}/peaks/all_peaks.tmp.bed
bedtools merge -i ${DATA_PATH}/peaks/all_peaks.tmp.bed > ${DATA_PATH}/peaks/all_peaks.tmp.merged.bed
cat ${DATA_PATH}/peaks/all_peaks.tmp.merged.bed | awk '{ print $0"\tpeak_"NR }' > ${DATA_PATH}/peaks/all_peaks.canonical.replicated.no_blacklist.bed
rm ${DATA_PATH}/peaks/all_peaks.tmp.bed ${DATA_PATH}/peaks/all_peaks.tmp.merged.bed

# Lift over merged peak set
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap1.chain \
		${DATA_PATH}/peaks/all_peaks.canonical.replicated.no_blacklist.bed \
		${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap1_lifted.bed
	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap2.chain \
		${DATA_PATH}/peaks/all_peaks.canonical.replicated.no_blacklist.bed \
		${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap2_lifted.bed

	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap1_lifted.bed |\
		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | sort -k1,1 -k2,2n > ${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap1.bed
	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap2_lifted.bed |\
		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | sort -k1,1 -k2,2n > ${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap2.bed

	rm ${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap*_lifted.bed
	rm ${DATA_PATH}/peaks/all_peaks_${SAMPLE}_hap*_lifted.bed.unmap
done

# Lift over individual peak sets
#for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
#do
#	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap1.chain \
#		${DATA_PATH}/peaks/${SAMPLE}_peaks.canonical.replicated.no_blacklist.bed \
#		${DATA_PATH}/peaks/${SAMPLE}_peaks_hap1_lifted.bed
#	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap2.chain \
#		${DATA_PATH}/peaks/${SAMPLE}_peaks.canonical.replicated.no_blacklist.bed \
#		${DATA_PATH}/peaks/${SAMPLE}_peaks_hap2_lifted.bed
#
#	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/peaks/${SAMPLE}_peaks_hap1_lifted.bed |\
#		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | sort -k1,1 -k2,2n > ${DATA_PATH}/peaks/${SAMPLE}_peaks_hap1.bed
#	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/peaks/${SAMPLE}_peaks_hap2_lifted.bed |\
#		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | sort -k1,1 -k2,2n > ${DATA_PATH}/peaks/${SAMPLE}_peaks_hap2.bed
#
#	rm ${DATA_PATH}/peaks/${SAMPLE}_peaks_hap*_lifted.bed
#	rm ${DATA_PATH}/peaks/${SAMPLE}_peaks_hap*_lifted.bed.unmap
#done
