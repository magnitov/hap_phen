eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate chrombpnet_gpu

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	chrombpnet prep nonpeaks \
		-g ./hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
		-p ./train_peaks/${SAMPLE}.pval_0.01.no_as_FDR_0.1.narrowPeak \
		-c ./hg38/hg38.chrom.sizes \
		-fl ./folds/fold_0.json \
		-br ./hg38/hg38-blacklist.v2.bed.gz \
		-o ./nonpeaks/${SAMPLE} \
		-s 42
done
