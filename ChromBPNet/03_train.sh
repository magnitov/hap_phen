#!/usr/bin/env bash
#SBATCH -p oncode
#SBATCH -t 7-1:00:00
#SBATCH --cpus-per-gpu=32
#SBATCH --gres=gpu:1
#SBATCH --mem=50000mb
#SBATCH --output=../logs/03_train.out
#SBATCH --error=../logs/03_train.err
#SBATCH -w wallace

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate chrombpnet_gpu

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	chrombpnet pipeline \
		-ibam ./bam/${SAMPLE}_merged_ref_Aligned.sortedByCoord.mapq.nodups.bam \
		-d "ATAC" \
		-g ./hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
		-c ./hg38/hg38.chrom.sizes \
		-p ./train_peaks/${SAMPLE}.pval_0.01.no_as_FDR_0.1.narrowPeak \
		-n ./nonpeaks/${SAMPLE}_negatives.bed \
		-fl ./folds/fold_0.json \
		-b ./bias_model/${SAMPLE}/models/${SAMPLE}_bias.h5 \
		-o ./model/${SAMPLE}/ \
		-fp ${SAMPLE}
done
