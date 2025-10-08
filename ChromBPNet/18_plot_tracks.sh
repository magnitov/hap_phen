INPUT_PATH=/DATA/users/m.magnitov/hap_phen/chromBPNet

###
### Examples of observed and predicted ATAC-seq profiles
###
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/predictions_chrombpnet_CD40.ini \
	-o ${INPUT_PATH}/example_tracks/predictions_chrombpnet_CD40.pdf --dpi 300 \
	--region chr20:46,118,100-46,118,350 --width 14 --height 21 --fontSize 8 --trackLabelHAlign right

pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/predictions_chrombpnet_HTRA4.ini \
	-o ${INPUT_PATH}/example_tracks/predictions_chrombpnet_HTRA4.pdf --dpi 300 \
	--region chr8:38,977,390-38,977,640 --width 14 --height 21 --fontSize 8 --trackLabelHAlign right

###
### Examples of ChIP-seq profiles at predicted motif disruptions 
###
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/motif_AP1_example_NA12878_ChIPseq.ini \
	-o ${INPUT_PATH}/example_tracks/motif_AP1_example_NA12878_ChIPseq.pdf --dpi 300 \
	--region chr5:138,599,500-138,603,500 --width 8 --height 10 --fontSize 8 --trackLabelHAlign right
    
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/motif_ETS_example_NA12878_ChIPseq.ini \
	-o ${INPUT_PATH}/example_tracks/motif_ETS_example_NA12878_ChIPseq.pdf --dpi 300 \
	--region chr5:110,714,600-110,718,600 --width 8 --height 10 --fontSize 8 --trackLabelHAlign right

###
### Examples of ATAC-seq and TT-seq profiles at linked genes 
###
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/gene_TBC1D4_example_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/gene_TBC1D4_example_ATACseq_TTseq.pdf --dpi 300 \
	--region chr13:75,265,000-75,500,000 --width 14 --height 24 --fontSize 8 --trackLabelHAlign right
    
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/gene_SLFN5_example_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/gene_SLFN5_example_ATACseq_TTseq.pdf --dpi 300 \
	--region chr17:35,235,000-35,280,000 --width 14 --height 21 --fontSize 8 --trackLabelHAlign right

pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/gene_CAV2_example_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/gene_CAV2_example_ATACseq_TTseq.pdf --dpi 300 \
	--region chr7:116,497,000-116,511,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right
    
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/gene_PAX8-AS1_example_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/gene_PAX8-AS1_example_ATACseq_TTseq.pdf --dpi 300 \
	--region chr2:113,233,000-113,278,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right
    
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/gene_PIK3R5_example_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/gene_PIK3R5_example_ATACseq_TTseq.pdf --dpi 300 \
	--region chr17:8,870,000-8,975,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right

###
### Examples of ATAC-seq and TT-seq profiles at balanced genes 
###
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/balanced_region_PPIF_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/balanced_region_PPIF_ATACseq_TTseq.pdf --dpi 300 \
	--region chr10:79,346,000-79,356,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right

pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/balanced_region_SDCCAG8_ATACseq_TTseq.ini \
	-o ${INPUT_PATH}/example_tracks/balanced_region_SDCCAG8_ATACseq_TTseq.pdf --dpi 300 \
	--region chr1:243,276,000-243,301,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right
