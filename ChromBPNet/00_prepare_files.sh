DATA_PATH=/DATA/users/m.magnitov/hap_phen

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	# Softlink genomes
	ln -s ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap1.fa  ${DATA_PATH}/chromBPNet/genomes/${SAMPLE}_hap1.fa
	ln -s ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap2.fa  ${DATA_PATH}/chromBPNet/genomes/${SAMPLE}_hap2.fa
	ln -s ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap1.chrom.sizes ${DATA_PATH}/chromBPNet/genomes/${SAMPLE}_hap1.chrom.sizes
	ln -s ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap2.chrom.sizes ${DATA_PATH}/chromBPNet/genomes/${SAMPLE}_hap2.chrom.sizes

	# Softlink BAM files
	ln -s ${DATA_PATH}/ATACseq/bam/${SAMPLE}_merged_ref_Aligned.sortedByCoord.mapq.nodups.bam ${DATA_PATH}/chromBPNet/bam/${SAMPLE}_merged_ref_Aligned.sortedByCoord.mapq.nodups.bam
done

ln -s /DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${DATA_PATH}/chromBPNet/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
ln -s /DATA/users/m.magnitov/genomes/hg38.chrom.sizes ${DATA_PATH}/chromBPNet/hg38/hg38.chrom.sizes
ln -s /DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed.gz ${DATA_PATH}/chromBPNet/hg38/hg38-blacklist.v2.bed.gz