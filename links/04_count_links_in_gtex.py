import os
import pandas as pd
import numpy as np
import bioframe

def read_variants(sample):
    variants = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/{sample}.phased.het.bed', 
                           sep = '\t', header = None)
    variants[11] = [int(x.split('|')[0]) for x in variants[10].values]
    variants[12] = [int(x.split('|')[1].split(':')[0]) for x in variants[10].values]
    variants = variants[[0, 2, 2, 3, 5, 6, 11, 12]]
    variants.columns = ['chrom', 'start', 'end', 'variant_id', 'ref', 'alt', 'allele_ref', 'allele_alt']
    return(variants)

def get_variants_for_links(phased_variants, links):
    phased_variants_overlap_links = bioframe.overlap(phased_variants, links,  
                                                     cols2 = ('chrom_peak', 'start_peak', 'end_peak'),
                                                     suffixes = ('_variant', ''))
    phased_variants_overlap_links = phased_variants_overlap_links.dropna()
    phased_variants_overlap_links = phased_variants_overlap_links.drop('distance', axis = 1)

    putative_non_overlaps = links[~links['peak_id'].isin(phased_variants_overlap_links['peak_id'])]['peak_id'].values
    if len(putative_non_overlaps) != 0:
        links_non_overlaps = links[links['peak_id'].isin(putative_non_overlaps)]

        phased_variants_overlap_non_overlaps = bioframe.closest(links_non_overlaps, phased_variants, 
                     cols1 = ('chrom_peak', 'start_peak', 'end_peak'),
                     suffixes = ('', '_variant'))
        sort_order = list(phased_variants_overlap_non_overlaps.columns[-9:]) + list(phased_variants_overlap_non_overlaps.columns[:-9])
        phased_variants_overlap_non_overlaps = phased_variants_overlap_non_overlaps[sort_order]
        phased_variants_overlap_non_overlaps = phased_variants_overlap_non_overlaps.drop('distance', axis = 1)

    phased_variants_overlap_links = pd.concat([phased_variants_overlap_links, phased_variants_overlap_non_overlaps], ignore_index = True)
    phased_variants_overlap_links = phased_variants_overlap_links.sort_values(['chrom_variant', 'start_variant'])
    phased_variants_overlap_links.index = np.arange(len(phased_variants_overlap_links))
    phased_variants_overlap_links['gene_id'] = [x.split('.')[0] for x in phased_variants_overlap_links['gene_id']]
    phased_variants_overlap_links['variant_link_id'] = ['variant_link_' + str(x) for x in np.arange(len(phased_variants_overlap_links))+1]
    return(phased_variants_overlap_links)

gtex_data_all = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/links/gtex_data/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Cells_EBV-transformed_lymphocytes.allpairs.txt.gz', sep = '\t')
gtex_data_all['chrom'] = [x.split('_')[0] for x in gtex_data_all['variant_id']]
gtex_data_all['start'] = [int(x.split('_')[1]) for x in gtex_data_all['variant_id']]
gtex_data_all['end'] = [int(x.split('_')[1]) for x in gtex_data_all['variant_id']]
gtex_data_all['ref'] = [x.split('_')[2] for x in gtex_data_all['variant_id']]
gtex_data_all['alt'] = [x.split('_')[3] for x in gtex_data_all['variant_id']]
gtex_data_all['gene_id'] = [x.split('.')[0] for x in gtex_data_all['gene_id']]
gtex_data_all = gtex_data_all[['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'maf', 'slope']]

gtex_data_sign = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/links/gtex_data/Cells_EBV-transformed_lymphocytes.v8.signif_variant_gene_pairs.txt.gz', sep = '\t')
gtex_data_sign['chrom'] = [x.split('_')[0] for x in gtex_data_sign['variant_id']]
gtex_data_sign['start'] = [int(x.split('_')[1]) for x in gtex_data_sign['variant_id']]
gtex_data_sign['end'] = [int(x.split('_')[1]) for x in gtex_data_sign['variant_id']]
gtex_data_sign['ref'] = [x.split('_')[2] for x in gtex_data_sign['variant_id']]
gtex_data_sign['alt'] = [x.split('_')[3] for x in gtex_data_sign['variant_id']]
gtex_data_sign['gene_id'] = [x.split('.')[0] for x in gtex_data_sign['gene_id']]
gtex_data_sign = gtex_data_sign[['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'maf', 'slope']]

gtex_data_variants = pd.read_csv('/DATA/users/m.magnitov/hap_phen/links/gtex_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz', sep = '\t')
gtex_data_variants['start'] = gtex_data_variants['variant_pos']
gtex_data_variants['end'] = gtex_data_variants['variant_pos']
gtex_data_variants = gtex_data_variants[['chr', 'start', 'end', 'ref', 'alt']]


gtex_abundance = pd.DataFrame()
for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:
    print(sample)
    phased_variants = read_variants(sample)

    links_allele_specific = pd.read_csv(f'./links/links_{sample}_allele_specific.txt', sep = '\t', header = 0)
    links_balanced = pd.read_csv(f'./links/links_{sample}_balanced.txt', sep = '\t', header = 0)

    variants_links_allele_specific = get_variants_for_links(phased_variants, links_allele_specific)
    variants_links_balanced = get_variants_for_links(phased_variants, links_balanced)

    variants_links_allele_specific_in_gtex_sign = variants_links_allele_specific.merge(gtex_data_sign, 
                                                                                       left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant', 'gene_id'], 
                                                                                       right_on = ['chrom', 'start', 'ref', 'alt', 'gene_id'])

    variants_links_allele_specific_in_gtex_all = variants_links_allele_specific.merge(gtex_data_all, 
                                                                                      left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant', 'gene_id'], 
                                                                                      right_on = ['chrom', 'start', 'ref', 'alt', 'gene_id'])

    variants_links_allele_specific_in_gtex_variants = variants_links_allele_specific.merge(gtex_data_variants, 
                                                                                           left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant'], 
                                                                                           right_on = ['chr', 'start', 'ref', 'alt'])

    variants_links_balanced_in_gtex_sign = variants_links_balanced.merge(gtex_data_sign, 
                                                                         left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant', 'gene_id'], 
                                                                         right_on = ['chrom', 'start', 'ref', 'alt', 'gene_id'])

    variants_links_balanced_in_gtex_all = variants_links_balanced.merge(gtex_data_all, 
                                                                        left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant', 'gene_id'], 
                                                                        right_on = ['chrom', 'start', 'ref', 'alt', 'gene_id'])

    variants_links_balanced_in_gtex_variants = variants_links_balanced.merge(gtex_data_variants, 
                                                                             left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant'], 
                                                                             right_on = ['chr', 'start', 'ref', 'alt'])
    
    gtex_sign_count, gtex_all_count, gtex_variants_count, other_count = 0, 0, 0, 0
    for variant_link_id in variants_links_allele_specific['variant_link_id'].values:
        if variant_link_id in variants_links_allele_specific_in_gtex_sign['variant_link_id'].values:
            gtex_sign_count += 1
        elif variant_link_id in variants_links_allele_specific_in_gtex_all['variant_link_id'].values:
            gtex_all_count += 1
        elif variant_link_id in variants_links_allele_specific_in_gtex_variants['variant_link_id'].values:
            gtex_variants_count += 1
        else:
            other_count += 1
    counts_total_as = np.array([gtex_sign_count, gtex_all_count, gtex_variants_count, other_count])
    print(counts_total_as)
    
    gtex_sign_count, gtex_all_count, gtex_variants_count, other_count = 0, 0, 0, 0
    for variant_link_id in variants_links_balanced['variant_link_id'].values:
        if variant_link_id in variants_links_balanced_in_gtex_sign['variant_link_id'].values:
            gtex_sign_count += 1
        elif variant_link_id in variants_links_balanced_in_gtex_all['variant_link_id'].values:
            gtex_all_count += 1
        elif variant_link_id in variants_links_balanced_in_gtex_variants['variant_link_id'].values:
            gtex_variants_count += 1
        else:
            other_count += 1
    counts_total_balanced = np.array([gtex_sign_count, gtex_all_count, gtex_variants_count, other_count])
    print(counts_total_balanced)

    gtex_abundance[sample + '_as'] = counts_total_as
    gtex_abundance[sample + '_balanced'] = counts_total_balanced

gtex_abundance.index = ['eqtl', 'tested_not_significant', 'not_tested_low_maf', 'missed']
gtex_abundance.to_csv('/DATA/users/m.magnitov/hap_phen/links/gtex_data/links_in_gtex_counts.txt', sep = '\t', index = 1, header = 1)
