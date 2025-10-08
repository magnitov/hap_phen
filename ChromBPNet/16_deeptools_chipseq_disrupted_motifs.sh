DATA_PATH=/DATA/users/m.magnitov/hap_phen/

mkdir -p ${DATA_PATH}/chromBPNet/deeptools_piles

#as balanced control RdPu BuGn Greys

###
### AP1: BATF, JUNB
###
for REGION_TYPE in control # as balanced control
do
	computeMatrix reference-point -p 32 --referencePoint center -a 2000 -b 2000 --missingDataAsZero --skipZeros \
		-S ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap1.hg38.bw ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap2.hg38.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_BATF_ENCSR000BGT_hap1.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_BATF_ENCSR000BGT_hap2.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_BATF_ENCSR000BGT.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_JUNB_ENCSR897MMC_hap1.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_JUNB_ENCSR897MMC_hap2.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_JUNB_ENCSR897MMC.hg38.input.bw \
		-R ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/AP1_hap1_${REGION_TYPE}.bed ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/AP1_hap2_${REGION_TYPE}.bed \
		-o ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_AP1_hits_${REGION_TYPE}.gz

	plotHeatmap -m ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_AP1_hits_${REGION_TYPE}.gz -out ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_AP1_hits_${REGION_TYPE}.pdf \
		--zMin 1 --zMax 200 200 14 14 14 8 8 8 \
		--yMin 0 --yMax 200 200 14 14 14 8 8 8 \
		--heatmapHeight 10 --heatmapWidth 3 --refPointLabel "" --xAxisLabel "" --legendLocation none --sortRegions keep \
		--samplesLabel ATAC_hap1 ATAC_hap2 BATF_hap1 BATF_hap2 BATF_full JUNB_hap1 JUNB_hap2 JUNB_full \
		--colorMap Blues Blues Greys Greys Greys Greys Greys Greys \
		--regionsLabel "${REGION_TYPE}_1" "${REGION_TYPE}_2"
done

###
### NFKB: RELA, RELB
###
for REGION_TYPE in control # as balanced control
do
	computeMatrix reference-point -p 32 --referencePoint center -a 2000 -b 2000 --missingDataAsZero --skipZeros \
		-S ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap1.hg38.bw ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap2.hg38.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_RELA_Zhao_hap1.hg38.rpgc.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_RELA_Zhao_hap2.hg38.rpgc.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_RELA_Zhao.hg38.rpgc.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_RELB_Zhao_hap1.hg38.rpgc.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_RELB_Zhao_hap2.hg38.rpgc.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_RELB_Zhao.hg38.rpgc.bw \
		-R ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/NFKB_hap1_${REGION_TYPE}.bed ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/NFKB_hap2_${REGION_TYPE}.bed \
		-o ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_NFKB_hits_${REGION_TYPE}.gz

	plotHeatmap -m ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_NFKB_hits_${REGION_TYPE}.gz -out ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_NFKB_hits_${REGION_TYPE}.pdf \
		--zMin 1 --zMax 200 200 450 450 20 750 750 40 \
		--yMin 0 --yMax 200 200 450 450 20 750 750 40 \
		--heatmapHeight 10 --heatmapWidth 3 --refPointLabel "" --xAxisLabel "" --legendLocation none --sortRegions keep \
		--samplesLabel ATAC_hap1 ATAC_hap2 RELA_hap1 RELA_hap2 RELA_full RELB_hap1 RELB_hap2 RELB_full \
		--colorMap Blues Blues Greys Greys Greys Greys Greys Greys \
		--regionsLabel "${REGION_TYPE}_1" "${REGION_TYPE}_2"
done

###
### RUNX: RUNX3
###
for REGION_TYPE in control # as balanced control
do
	computeMatrix reference-point -p 32 --referencePoint center -a 2000 -b 2000 --missingDataAsZero --skipZeros \
		-S ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap1.hg38.bw ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap2.hg38.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_RUNX3_ENCSR000BRI_hap1.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_RUNX3_ENCSR000BRI_hap2.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_RUNX3_ENCSR000BRI.hg38.input.bw \
		-R ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/RUNX_hap1_${REGION_TYPE}.bed ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/RUNX_hap2_${REGION_TYPE}.bed \
		-o ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_RUNX_hits_${REGION_TYPE}.gz

	plotHeatmap -m ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_RUNX_hits_${REGION_TYPE}.gz -out ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_RUNX_hits_${REGION_TYPE}.pdf \
		--zMin 1 --zMax 200 200 8 8 8 \
		--yMin 0 --yMax 200 200 8 8 8 \
		--heatmapHeight 10 --heatmapWidth 3 --refPointLabel "" --xAxisLabel "" --legendLocation none --sortRegions keep \
		--samplesLabel ATAC_hap1 ATAC_hap2 RUNX3_hap1 RUNX3_hap2 RUNX3_full \
		--colorMap Blues Blues Greys Greys Greys \
		--regionsLabel "${REGION_TYPE}_1" "${REGION_TYPE}_2"
done

###
### IRF: IRF4
###
for REGION_TYPE in control # as balanced control
do
	computeMatrix reference-point -p 32 --referencePoint center -a 2000 -b 2000 --missingDataAsZero --skipZeros \
		-S ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap1.hg38.bw ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap2.hg38.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_IRF4_ENCSR000BGY_hap1.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_IRF4_ENCSR000BGY_hap2.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_IRF4_ENCSR000BGY.hg38.input.bw \
		-R ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/IRF_hap1_${REGION_TYPE}.bed ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/IRF_hap2_${REGION_TYPE}.bed \
		-o ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_IRF_hits_${REGION_TYPE}.gz

	plotHeatmap -m ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_IRF_hits_${REGION_TYPE}.gz -out ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_IRF_hits_${REGION_TYPE}.pdf \
		--zMin 1 --zMax 200 200 5 5 5 \
		--yMin 0 --yMax 200 200 5 5 5 \
		--heatmapHeight 10 --heatmapWidth 3 --refPointLabel "" --xAxisLabel "" --legendLocation none --sortRegions keep \
		--samplesLabel ATAC_hap1 ATAC_hap2 IRF4_hap1 IRF4_hap2 IRF4_full \
		--colorMap Blues Blues Greys Greys Greys \
		--regionsLabel "${REGION_TYPE}_1" "${REGION_TYPE}_2"
done

###
### ETS: SPI1
###
for REGION_TYPE in control # as balanced control
do
	computeMatrix reference-point -p 32 --referencePoint center -a 2000 -b 2000 --missingDataAsZero --skipZeros \
		-S ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap1.hg38.bw ${DATA_PATH}/ATACseq/bigwig_haplotypes/ATACseq_NA12878_hap2.hg38.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_SPI1_ENCSR000BGQ_hap1.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig_haplotypes/ChIPseq_SPI1_ENCSR000BGQ_hap2.hg38.input.bw ${DATA_PATH}/ChIPseq/bigwig/ChIPseq_SPI1_ENCSR000BGQ.hg38.input.bw \
		-R ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/ETS_hap1_${REGION_TYPE}.bed ${DATA_PATH}/chromBPNet/variants_annotation/motif_hits_NA12878/ETS_hap2_${REGION_TYPE}.bed \
		-o ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_ETS_hits_${REGION_TYPE}.gz

	plotHeatmap -m ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_ETS_hits_${REGION_TYPE}.gz -out ${DATA_PATH}/chromBPNet/deeptools_piles/ChIPseq_ETS_hits_${REGION_TYPE}.pdf \
		--zMin 1 --zMax 200 200 10 10 10 \
		--yMin 0 --yMax 200 200 10 10 10 \
		--heatmapHeight 10 --heatmapWidth 3 --refPointLabel "" --xAxisLabel "" --legendLocation none --sortRegions keep \
		--samplesLabel ATAC_hap1 ATAC_hap2 SPI1_hap1 SPI1_hap2 SPI1_full \
		--colorMap Blues Blues Greys Greys Greys \
		--regionsLabel "${REGION_TYPE}_1" "${REGION_TYPE}_2"
done