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

def get_chrom_size(chromosome):
    with open(f'/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0/fasta/genome.dict', 'r') as f:
        chromsizes = f.readlines()
    chromsizes = chromsizes[1:]
    chromnames = [x.split('\t')[1].split(':')[1] for x in chromsizes]
    chromsizes = [int(x.split('\t')[2].split(':')[1]) for x in chromsizes]
    
    chromsizes = pd.DataFrame({'chr': chromnames, 'size': chromsizes})
    
    return(chromsizes[chromsizes['chr'] == 'chr' + chromosome]['size'].values[0])

def get_haploblocks_data(sample, chromosome):
    with open(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/output/chr{chromosome}.hap', 'r') as f:
        haplotypes_stats = f.readlines()
    haplotypes_stats = [x for x in haplotypes_stats if 'BLOCK' in x]

    snp_num, snp_phased_num, block_span = [], [], []
    for block in haplotypes_stats:
        snp_num.append(int(block.split('len:')[1].split('phased:')[0][1:-1]))
        snp_phased_num.append(int(block.split('phased:')[1].split('SPAN:')[0][1:-1]))
        block_span.append(int(block.split('SPAN:')[1].split('fragment')[0][1:-1]))
        
    haploblocks_data = pd.DataFrame([snp_phased_num, block_span]).transpose()
    haploblocks_data.columns = ['het_snp_phased', 'span']
    haploblocks_data = haploblocks_data.sort_values(['span'], ascending = 0)
    
    chrom_size = get_chrom_size(chromosome)
    het_variants = get_number_of_hetvariants_per_chrom_from_hapcut2(sample, chromosome)
    
    haploblocks_data['het_snp_phased_frac'] = [x/het_variants for x in haploblocks_data['het_snp_phased']]
    haploblocks_data['span_frac'] = [x/chrom_size for x in haploblocks_data['span']]
    
    haploblocks_data['mvp'] = haploblocks_data['span_frac']*haploblocks_data['het_snp_phased_frac']
    
    return(haploblocks_data.sort_values(['mvp'], ascending = False))

def get_index_of_mvp_in_vcf(sample, chromosome, blocks):
    
    num_snps_mvp = blocks['het_snp_phased'].values[0]

    vcf = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/output/chr{chromosome}.hap.phased.VCF', 
                      comment = '#', delim_whitespace = True, header = None)
    
    snps_per_block = pd.DataFrame(np.unique([x.split(':')[-1] for x in vcf[9]], return_counts = 1)).transpose()
    snps_per_block[2] = abs(snps_per_block[1] - num_snps_mvp)
    snps_per_block = snps_per_block.sort_values([2], ascending = True)
    
    return(snps_per_block[0].values[0])

def create_mvp_vcf(sample, chromosome, block_index):
    data_path = '/DATA/users/m.magnitov/hap_phen/personal_genomes'
    
    os.system('grep' + " '^#' " +  data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.VCF' +\
              ' > ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp_header.VCF')
    os.system('grep :' + block_index + '$ ' +  data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.VCF' +\
              ' | grep -e' + " '1|0' -e" + " '0|1'" +\
              ' > ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp_snps.VCF')
    os.system('cat ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp_header.VCF' +\
              ' ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp_snps.VCF' +\
              ' > ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp.VCF')
    os.system('rm ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp_header.VCF' +\
              ' ' + data_path + '/' + sample + '/output/chr' + chromosome + '.hap.phased.mvp_snps.VCF')
    
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
        blocks = get_haploblocks_data(sample, chromosome)
        block_index_vcf = get_index_of_mvp_in_vcf(sample, chromosome, blocks)
        create_mvp_vcf(sample, chromosome, block_index_vcf)
    # Chromosome X
    if sample in ['NA12878', 'NA19238', 'NA19240', 'HG00513', 'HG00514', 
                  'HG00732', 'HG00733', 'HG03486', 'HG02601', 'HG03464']:
        blocks = get_haploblocks_data(sample, 'X')
        block_index_vcf = get_index_of_mvp_in_vcf(sample, 'X', blocks)
        create_mvp_vcf(sample, 'X', block_index_vcf)
