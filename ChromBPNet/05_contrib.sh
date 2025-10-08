#!/usr/bin/env bash
#SBATCH -p oncode
#SBATCH -t 7-1:00:00
#SBATCH --cpus-per-gpu=32
#SBATCH --gres=gpu:1
#SBATCH --mem=50000mb
#SBATCH --output=../logs/05_contrib.out
#SBATCH --error=../logs/05_contrib.err
#SBATCH -w wallace

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate chrombpnet_gpu

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	chrombpnet contribs_bw \
		-m ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet_nobias.h5 \
		-g ./genomes/${SAMPLE}_hap1.fa \
		-c ./genomes/${SAMPLE}_hap1.chrom.sizes \
		-r ./predict_peaks/${SAMPLE}_hap1.narrowPeak \
		-pc profile \
		-op ./contrib/${SAMPLE}_hap1

	chrombpnet contribs_bw \
		-m ./model/${SAMPLE}/models/${SAMPLE}_chrombpnet_nobias.h5 \
		-g ./genomes/${SAMPLE}_hap2.fa \
		-c ./genomes/${SAMPLE}_hap2.chrom.sizes \
		-r ./predict_peaks/${SAMPLE}_hap2.narrowPeak \
		-pc profile \
		-op ./contrib/${SAMPLE}_hap2
done
