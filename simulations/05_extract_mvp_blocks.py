import os
import pandas as pd
import numpy as np

def get_haploblocks_data(lr_type, depth_hic, depth_lr, sampling_num):
    hap_path = '/DATA/users/m.magnitov/hap_phen/simulations/haplotypes'
    with open(f'{hap_path}/output_hic_{lr_type}/chr1_{depth_hic}x_hic_{depth_lr}x_{lr_type}_ds{sampling_num}.hap', 'r') as f:
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
    
    # size of chr1 in hg38 and number if het SNPs for chr1 in NA12878
    chrom_size = 248956422
    het_variants = 195087
    
    haploblocks_data['het_snp_phased_frac'] = [x/het_variants for x in haploblocks_data['het_snp_phased']]
    haploblocks_data['span_frac'] = [x/chrom_size for x in haploblocks_data['span']]
    
    # most variants phased block is defined by maximizing span annd phased het SNPs
    haploblocks_data['mvp'] = haploblocks_data['span_frac']*haploblocks_data['het_snp_phased_frac']
    
    return(haploblocks_data.sort_values(['mvp'], ascending = False))

def get_index_of_mvp_in_vcf(lr_type, depth_hic, depth_lr, sampling_num, blocks): 
    hap_path = '/DATA/users/m.magnitov/hap_phen/simulations/haplotypes'
    vcf = pd.read_csv(f'{hap_path}/output_hic_{lr_type}/chr1_{depth_hic}x_hic_{depth_lr}x_{lr_type}_ds{sampling_num}.hap.phased.VCF', 
                      comment = '#', delim_whitespace = True, header = None)
    vcf = vcf[vcf[9].str.contains('/') == False] 
    
    num_snps_mvp = blocks['het_snp_phased'].values[0]

    snps_per_block = pd.DataFrame(np.unique([x.split(':')[-1] for x in vcf[9]], return_counts = 1)).transpose()
    snps_per_block[2] = abs(snps_per_block[1] - num_snps_mvp)
    snps_per_block = snps_per_block.sort_values([2], ascending = True)
    
    return(snps_per_block[0].values[0])

def create_mvp_vcf(lr_type, depth_hic, depth_lr, sampling_num, block_index):
    hap_path = '/DATA/users/m.magnitov/hap_phen/simulations/haplotypes'
    
    os.system('grep' + " '^#' " +  hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.VCF' +\
              ' > ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp_header.VCF')
    os.system('grep :' + block_index + '$ ' +  hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.VCF' +\
              ' > ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp_snps.VCF')
    os.system('cat ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp_header.VCF' +\
              ' ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp_snps.VCF' +\
              ' > ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp.VCF')
    os.system('rm ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp_header.VCF' +\
              ' ' + hap_path + '/output_hic_' + lr_type +\
              '/chr1_' + depth_hic + 'x_hic_' + depth_lr + 'x_' + lr_type + '_ds' + sampling_num + '.hap.phased.mvp_snps.VCF')
    
    return(True)

### Create VCF files for MVP blocks
for lr_type in ['10X', 'tellseq', 'stlfr']:
    for sampling_num in np.arange(1, 21, 1):
        for depth_hic in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 32]:
            for depth_lr in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 32]:
            
                blocks = get_haploblocks_data(lr_type, depth_hic, depth_lr, sampling_num)
                block_index_vcf = get_index_of_mvp_in_vcf(lr_type, depth_hic, depth_lr, sampling_num, blocks)
                create_mvp_vcf(lr_type, str(depth_hic), str(depth_lr), str(sampling_num), block_index_vcf)
