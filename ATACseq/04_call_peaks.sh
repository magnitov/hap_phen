# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
BLACKLIST=/DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed

mkdir -p ${DATA_PATH}/peaks

# Call peaks
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	# Replicated
	for REPLICATE in rep1 rep2 rep3
	do
		macs2 callpeak -f BAMPE --gsize=hs --pvalue=0.01 --nomodel --keep-dup all \
			-t ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.nodups.bam \
			--name=${SAMPLE}_${REPLICATE} --outdir=${DATA_PATH}/peaks \
			2> ${DATA_PATH}/peaks/macs_${SAMPLE}_${REPLICATE}.log
	done

	# Pooled
	macs2 callpeak -f BAMPE --gsize=hs --pvalue=0.01 --nomodel --keep-dup all \
		-t ${DATA_PATH}/bam/${SAMPLE}_*_ref_Aligned.sortedByCoord.mapq.nodups.bam \
		--name=${SAMPLE}_merged --outdir=${DATA_PATH}/peaks \
		2> ${DATA_PATH}/peaks/macs_${SAMPLE}_merged.log
done

# Get replicated peaks and remove peaks from blacklisted regions
cd ${DATA_PATH}/peaks
rm *.xls *_summits.bed

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	# Remove non-canonical chromosomes
	grep -v -e 'chrUn' -e 'random' -e 'chrY' -e 'chrM' -e 'chrEBV' ${SAMPLE}_rep1_peaks.narrowPeak > ${SAMPLE}_rep1_peaks.canonical.bed
	grep -v -e 'chrUn' -e 'random' -e 'chrY' -e 'chrM' -e 'chrEBV' ${SAMPLE}_rep2_peaks.narrowPeak > ${SAMPLE}_rep2_peaks.canonical.bed
	grep -v -e 'chrUn' -e 'random' -e 'chrY' -e 'chrM' -e 'chrEBV' ${SAMPLE}_rep3_peaks.narrowPeak > ${SAMPLE}_rep3_peaks.canonical.bed
	grep -v -e 'chrUn' -e 'random' -e 'chrY' -e 'chrM' -e 'chrEBV' ${SAMPLE}_merged_peaks.narrowPeak > ${SAMPLE}_merged_peaks.canonical.bed
	
	# Get replicated peaks set
	bedtools intersect -wo -a ${SAMPLE}_merged_peaks.canonical.bed -b ${SAMPLE}_rep1_peaks.canonical.bed | awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-4 | sort -k1,1 -k2,2n | uniq > overlap_rep1.bed
	bedtools intersect -wo -a ${SAMPLE}_merged_peaks.canonical.bed -b ${SAMPLE}_rep2_peaks.canonical.bed | awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-4 | sort -k1,1 -k2,2n | uniq > overlap_rep2.bed              
	bedtools intersect -wo -a ${SAMPLE}_merged_peaks.canonical.bed -b ${SAMPLE}_rep3_peaks.canonical.bed | awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-4 | sort -k1,1 -k2,2n | uniq > overlap_rep3.bed
	cat overlap_rep1.bed overlap_rep2.bed overlap_rep3.bed | sort -k1,1 -k2,2n | uniq -c | awk '{ if ($1>=2) print $2"\t"$3"\t"$4 }' > ${SAMPLE}_peaks.canonical.replicated.bed

	# Remove peaks from blacklisted regions
	bedtools intersect -v -a ${SAMPLE}_peaks.canonical.replicated.bed -b ${BLACKLIST} | awk -v name="$SAMPLE" '{ print $0"\t"name"_peak_"NR }' > ${SAMPLE}_peaks.canonical.replicated.no_blacklist.bed
	
	# Remove intermediate files	
	rm overlap_rep*.bed ${SAMPLE}_*_peaks.canonical.bed ${SAMPLE}_peaks.canonical.replicated.bed
done

# Merge all peaks from five cell lines into unified peak set
cat ${DATA_PATH}/peaks/*_peaks.canonical.replicated.no_blacklist.bed | sort -k1,1 -k2,2n > ${DATA_PATH}/peaks/all_peaks.tmp.bed
bedtools merge -i ${DATA_PATH}/peaks/all_peaks.tmp.bed > ${DATA_PATH}/peaks/all_peaks.canonical.replicated.no_blacklist.bed
rm ${DATA_PATH}/peaks/all_peaks.tmp.bed