HAPLOTYPES_PATH=/DATA/users/m.magnitov/hap_phen/simulations/haplotypes
PLATINUM_VCF=/DATA/users/m.magnitov/hap_phen/variants/NA12878/NA12878.platinum.phased.chr1.het.vcf
CALC_HAP_STATS=/DATA/users/m.magnitov/software/HapCUT2-1.3.3/utilities/calculate_haplotype_statistics.py

for SAMPLING_NUM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	echo $SAMPLING_NUM
	for DEPTH_HiC in 1 2 3 4 5 6 7 8 10 12 16 20 24 32
	do
		for DEPTH_LR in 1 2 3 4 5 6 7 8 10 12 16 20 24 32
		do
			python ${CALC_HAP_STATS} -i -v1 ${HAPLOTYPES_PATH}/output_hic_10X/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_10X_ds${SAMPLING_NUM}.hap.phased.mvp.VCF -v2 ${PLATINUM_VCF} > ${HAPLOTYPES_PATH}/output_hic_10X/hap_stats_${DEPTH_HiC}x_hic_${DEPTH_LR}x_10X_ds${SAMPLING_NUM}.txt
			python ${CALC_HAP_STATS} -i -v1 ${HAPLOTYPES_PATH}/output_hic_tellseq/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_tellseq_ds${SAMPLING_NUM}.hap.phased.mvp.VCF -v2 ${PLATINUM_VCF} > ${HAPLOTYPES_PATH}/output_hic_tellseq/hap_stats_${DEPTH_HiC}x_hic_${DEPTH_LR}x_tellseq_ds${SAMPLING_NUM}.txt
			python ${CALC_HAP_STATS} -i -v1 ${HAPLOTYPES_PATH}/output_hic_stlfr/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_stlfr_ds${SAMPLING_NUM}.hap.phased.mvp.VCF -v2 ${PLATINUM_VCF} > ${HAPLOTYPES_PATH}/output_hic_stlfr/hap_stats_${DEPTH_HiC}x_hic_${DEPTH_LR}x_stlfr_ds${SAMPLING_NUM}.txt
		done
	done
done

