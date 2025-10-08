import os
import pandas as pd
import numpy as np

def get_number_of_hetvariants_per_chrom_from_hapcut2(sample, chromosome):
    with open(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/hapcut2_{sample}.log', 'r') as f:
        hapcut2_log = f.readlines()
    hapcut2_log = [x for x in hapcut2_log if 'vcffile /DATA/users/m.magnitov/hap_phen/variants/' + sample + '/split/chr' + chromosome + '.vcf' in x]
    hapcut2_log = np.unique(hapcut2_log)[0]
    hapcut2_log = int(hapcut2_log.split('variants')[2].lstrip().rstrip())
    return(hapcut2_log)

def get_chrom_size_for_assembly(assembly, chromosome):
    with open(f'/DATA/users/m.magnitov/genomes/refdata-{assembly}-2.1.0/fasta/genome.dict', 'r') as f:
        chromsizes = f.readlines()
    chromsizes = chromsizes[1:]
    chromnames = [x.split('\t')[1].split(':')[1] for x in chromsizes]
    chromsizes = [int(x.split('\t')[2].split(':')[1]) for x in chromsizes]
    
    chromsizes = pd.DataFrame({'chr': chromnames, 'size': chromsizes})
    
    return(chromsizes[chromsizes['chr'] == 'chr' + chromosome]['size'].values[0])

def get_haploblocks_data(sample, chromosome, assembly):
    phasing_10x = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/10X/{sample}/outs/phased_variants.vcf.gz', 
                              comment = '#', delim_whitespace = True, header = None, low_memory = False)
    phasing_10x['phase'] = [x.split(':')[0] for x in phasing_10x[9]]
    phasing_10x['index'] = [x.split(':')[2] for x in phasing_10x[9]]
    
    phasing_10x = phasing_10x[(phasing_10x['phase'] == '1|0') | (phasing_10x['phase'] == '0|1')]
    
    phasing_10x_chr = phasing_10x[phasing_10x[0] == 'chr' + chromosome]
    chrom_size = get_chrom_size_for_assembly(assembly, chromosome)
    het_variants = get_number_of_hetvariants_per_chrom_from_hapcut2(sample, chromosome)
    
    phased_hetvariants, phased_spans, indexes = [], [], []
    for block_index in np.unique(phasing_10x_chr['index']):
        phasing_10x_chr_index = phasing_10x_chr[phasing_10x_chr['index'] == block_index]
        
        indexes.append(block_index)
        phased_hetvariants.append(len(phasing_10x_chr_index))
        phased_spans.append(phasing_10x_chr_index[1].values[-1] - phasing_10x_chr_index[1].values[0])
        
    haploblocks_data = pd.DataFrame([indexes, phased_hetvariants, phased_spans]).transpose()
    haploblocks_data.columns = ['block_index', 'het_snp_phased', 'span']
    
    haploblocks_data['het_snp_phased_frac'] = [x/het_variants for x in haploblocks_data['het_snp_phased']]
    haploblocks_data['span_frac'] = [x/chrom_size for x in haploblocks_data['span']]
    
    haploblocks_data['mvp'] = haploblocks_data['span_frac']*haploblocks_data['het_snp_phased_frac']
    
    return(haploblocks_data.sort_values(['mvp'], ascending = False))

def create_mvp_vcf(sample, chromosome, block_index):
    data_path = '/DATA/users/m.magnitov/hap_phen/10X'
    
    os.system('zcat ' + data_path + '/' + sample + '/outs/phased_variants.vcf.gz | grep ' + " '^#' " +\
              ' > ' + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp_header.VCF')
    os.system('zcat ' + data_path + '/' + sample + '/outs/phased_variants.vcf.gz | grep :' + block_index + ': | ' +\
              'grep -e' + " '1|0' -e" + " '0|1' | " +\
              "awk '{ if ($1==" +'"chr' + chromosome + '"' + ") print $0 }' > " + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp_snps.VCF')
    os.system('cat ' + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp_header.VCF' +\
              ' ' + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp_snps.VCF' +\
              ' > ' + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp.VCF')
    os.system('rm ' + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp_header.VCF' +\
              ' ' + data_path + '/' + sample + '/chr' + chromosome + '.phased.mvp_snps.VCF')
    
    return(True)

### Create VCF files for MVP blocks
for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464',
               'NA19238', 'NA19239', 'NA19240', 
               'HG00512', 'HG00513', 'HG00514', 
               'HG00731', 'HG00732', 'HG00733',
               'NA20509', 'HG03486']:
    # Autosomes
    for chromosome in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
                       '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']:
        blocks = get_haploblocks_data(sample, chromosome, 'hg38')
        block_index_vcf = blocks['block_index'].values[0]
        create_mvp_vcf(sample, chromosome, block_index_vcf)
    # Chromosome X
    if sample in ['NA12878', 'NA12892', 'NA19238', 'NA19240', 'HG00513', 'HG00514', 
                  'HG00732', 'HG00733', 'HG03486', 'HG02601', 'HG03464']:
        blocks = get_haploblocks_data(sample, 'X', 'hg38')
        block_index_vcf = blocks['block_index'].values[0]
        create_mvp_vcf(sample, 'X', block_index_vcf)
