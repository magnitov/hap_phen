DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
CHIP_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq

mkdir -p ${DATA_PATH}/deeptools_piles

awk '{ if ($6<0) print $1"\t"$2"\t"$3 }' ${DATA_PATH}/asocr/NA12878_allele_specific.bed | grep -v 'chrX' > ${DATA_PATH}/deeptools_piles/NA12878_allele_specific_hap1.bed
awk '{ if ($6>0) print $1"\t"$2"\t"$3 }' ${DATA_PATH}/asocr/NA12878_allele_specific.bed | grep -v 'chrX' > ${DATA_PATH}/deeptools_piles/NA12878_allele_specific_hap2.bed

computeMatrix reference-point -p 32 --referencePoint center -a 3000 -b 3000 --missingDataAsZero --skipZeros \
	-S ${DATA_PATH}/bigwig_haplotypes/ATACseq_NA12878_hap1.hg38.bw ${DATA_PATH}/bigwig_haplotypes/ATACseq_NA12878_hap2.hg38.bw \
	${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K27ac_Kasowski_hap1.hg38.input.bw ${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K27ac_Kasowski_hap2.hg38.input.bw \
	${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K4me3_Kasowski_hap1.hg38.input.bw ${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K4me3_Kasowski_hap2.hg38.input.bw \
	${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K4me1_Kasowski_hap1.hg38.input.bw ${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K4me1_Kasowski_hap2.hg38.input.bw \
	${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K27me3_Kasowski_hap1.hg38.input.bw ${CHIP_PATH}/bigwig_haplotypes/ChIPseq_H3K27me3_Kasowski_hap2.hg38.input.bw \
	-R ${DATA_PATH}/deeptools_piles/NA12878_allele_specific_hap1.bed ${DATA_PATH}/deeptools_piles/NA12878_allele_specific_hap2.bed \
	-o ${DATA_PATH}/deeptools_piles/ChIPseq_signal_in_ATACseq_peaks.gz

plotHeatmap -m ${DATA_PATH}/deeptools_piles/ChIPseq_signal_in_ATACseq_peaks.gz -out ${DATA_PATH}/deeptools_piles/ChIPseq_signal_in_ATACseq_peaks.pdf \
	--zMin 0 --zMax 125 125 5 5 7 7 4 4 2 2 --yMin 0 --yMax 125 125 5 5 7 7 4 4 2 2 --heatmapHeight 12 --heatmapWidth 3 \
	--regionsLabel "Haplotype 1 (N=444)" "Haplotype 2 (N=439)" --refPointLabel "" --xAxisLabel "" --legendLocation none \
	--samplesLabel ATAC_hap1 ATAC_hap2 H3K27ac_hap1 H3K27ac_hap2 H3K4me3_hap1 H3K4me3_hap2 H3K4me1_hap1 H3K4me1_hap2 H3K27me3_hap1 H3K27me3_hap2 \
	--colorMap Blues Blues Oranges Oranges Reds Reds Purples Purples Greens Greens --sortUsingSamples 1 2
