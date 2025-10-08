#!/usr/bin/env bash
#SBATCH -p oncode
#SBATCH -t 3-1:00:00
#SBATCH --cpus-per-gpu=32
#SBATCH --gres=gpu:1
#SBATCH --mem=50000mb
#SBATCH --output=../logs/04_pred.out
#SBATCH --error=../logs/04_pred.err
#SBATCH -w wallace

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate chrombpnet_gpu

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	chrombpnet pred_bw \
		-cm ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet.h5
		-cmb ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet_nobias.h5 \
		-g ./genomes/${SAMPLE}_hap1.fa \
		-c ./genomes/${SAMPLE}_hap1.chrom.sizes \
		-r ./predict_peaks/${SAMPLE}_hap1.narrowPeak \
		-op ./pred/${SAMPLE}_hap1

	chrombpnet pred_bw \
		-cm ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet.h5
		-cmb ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet_nobias.h5 \
		-g ./genomes/${SAMPLE}_hap2.fa \
		-c ./genomes/${SAMPLE}_hap2.chrom.sizes \
		-r ./predict_peaks/${SAMPLE}_hap2.narrowPeak \
		-op ./pred/${SAMPLE}_hap2

	chrombpnet pred_bw \
		-cm ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet.h5
		-cmb ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet_nobias.h5 \
		-g ./hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
		-c ./hg38/hg38.chrom.sizes \
		-r ./predict_peaks/${SAMPLE}_hg38.narrowPeak \
		-op ./pred/${SAMPLE}_hg38
done

