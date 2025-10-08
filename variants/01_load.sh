#!/bin/bash

# individuals from 1000 Genomes
for CHROM in {1..22}
do
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
done
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
bcftools concat --threads 64 -O z -o 1kg_20220422.vcf.gz 1kGP_*.vcf.gz
rm 1kGP_*

# NA12878 from Platinum Genomes phasing
wget -O NA12878.platinum.phased.vcf.gz ftp://platgene_ro:''@ussd-ftp.illumina.com/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz

# LCLs from 1000 Genomes HGSVC phasing
wget -O HGSVC_20200814_freeze3.SNV.phased.vcf.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_snv_snv_alt.vcf.gz
wget -O HGSVC_20200814_freeze3.INDEL.phased.vcf.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_indel_insdel_alt.vcf.gz
