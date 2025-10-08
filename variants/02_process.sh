#!/bin/bash
VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants
GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

###
### preprocess vcfs
###
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	bcftools view --threads 64 --trim-alt-alleles -s ${SAMPLE} -O z -o ${SAMPLE}.vcf.gz 1kg_20220422.vcf.gz
	java -jar /DATA/users/m.magnitov/software/snpEff_v4_3p/SnpSift.jar filter -v "(countVariant() > 0 )" ${SAMPLE}.vcf.gz > ${SAMPLE}.filter.vcf
	cat ${SAMPLE}.filter.vcf | grep -v 'SVTYPE=' > ${SAMPLE}.filter.no_sv.vcf
	sed -i 's/Type=Integer,Number=1/Number=1,Type=Integer/g' ${SAMPLE}.filter.no_sv.vcf
	gzip ${SAMPLE}.filter.no_sv.vcf
	rm ${SAMPLE}.filter.vcf
	mv ${SAMPLE}.filter.no_sv.vcf.gz ${SAMPLE}.filter.vcf.gz
done

###
### rename chroms
###
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	bcftools annotate --threads 64 --rename-chrs chroms_to_rename.txt -O z -o ${SAMPLE}.filter.prefix.vcf.gz ${SAMPLE}.filter.vcf.gz
done

###
### uniform sites
###
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	###
	### normalize and left-align indels and to split the multiallelic sites into multiple rows
	###
	bcftools norm --threads 64 -m -any -f ${GENOME} -O z \
		-o ${SAMPLE}.filter.prefix.norm.vcf.gz \
		${SAMPLE}.filter.prefix.vcf.gz
	tabix -p vcf ${SAMPLE}.filter.prefix.norm.vcf.gz
    
	###
	### decompose complex variants and sort and unify resulting variants
	###
	vcfallelicprimitives ${SAMPLE}.filter.prefix.norm.vcf.gz --keep-info --keep-geno |\
		vt sort - | vt uniq - | bgzip -c > ${SAMPLE}.filter.prefix.norm.aprimitives.vcf.gz
	rm ${SAMPLE}.filter.prefix.norm.vcf.gz.tbi

	###
	### merge multiallelic sites into single rows
	###
	bcftools norm --threads 64 -m +any -f ${GENOME} -O z \
		-o ${SAMPLE}.filter.prefix.norm.aprimitives.merged.vcf.gz \
		${SAMPLE}.filter.prefix.norm.aprimitives.vcf.gz
   
	###
	### discard multiallelic sites
	###
	bcftools view --threads 64 -m2 -M2 -O v \
		-o ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.vcf \
		${SAMPLE}.filter.prefix.norm.aprimitives.merged.vcf.gz

	###
	### remove overlapping variants
	###
	vt remove_overlap -o ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_over_marked.vcf \
		${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.vcf
	java -jar /DATA/users/m.magnitov/software/snpEff_v4_3p/SnpSift.jar filter "(FILTER = 'PASS')" \
		${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_over_marked.vcf >\
		${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf
	gzip ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.vcf
	bgzip ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf
	rm ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_over_marked.vcf
done

###
### postprocess vcf
###
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	mkdir -p ${VCF_PATH}/${SAMPLE}
	mv ${VCF_PATH}/${SAMPLE}.* ${VCF_PATH}/${SAMPLE}/
	cd ${VCF_PATH}/${SAMPLE}

	###
	### split vcf by chromosome
	###
	tabix -p vcf ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf.gz
	mkdir -p ./split
	for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
	do
		bcftools view ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf.gz --regions chr${CHROM} -O v -o ./split/chr${CHROM}.vcf
	done
	rm ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf.gz.tbi
   
	###
	### softlinks for phasing
	###
	gunzip ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf.gz
	ln -s ${VCF_PATH}/${SAMPLE}/${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf ${VCF_PATH}/${SAMPLE}/${SAMPLE}.vcf

	###
	### autosomal variants for phasing statistics calculation
	###
	java -jar /DATA/users/m.magnitov/software/snpEff_v4_3p/SnpSift.jar filter -n "( CHROM = 'chrX' )" ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.vcf > ${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.autosomal.vcf
	ln -s ${VCF_PATH}/${SAMPLE}/${SAMPLE}.filter.prefix.norm.aprimitives.merged.biallelic.no_overlap.autosomal.vcf ${VCF_PATH}/${SAMPLE}/${SAMPLE}.autosomal.vcf
done

