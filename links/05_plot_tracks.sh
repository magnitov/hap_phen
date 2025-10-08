INPUT_PATH=/DATA/users/m.magnitov/hap_phen/links

###
### Examples of ATAC-seq and TT-seq profiles at linked genes 
###
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/link_ATACseq_TTseq_DAPP1.ini \
	-o ${INPUT_PATH}/example_tracks/link_ATACseq_TTseq_DAPP1.pdf --dpi 300 \
	--region chr4:99,809,000-99,878,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right
    
pyGenomeTracks --tracks ${INPUT_PATH}/example_tracks/link_ATACseq_TTseq_ZNF608.ini \
	-o ${INPUT_PATH}/example_tracks/link_ATACseq_TTseq_ZNF608.pdf --dpi 300 \
	--region chr5:124,615,000-124,760,000 --width 14 --height 15 --fontSize 8 --trackLabelHAlign right
