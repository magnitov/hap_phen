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

def read_significant_gtex_data(tissue_name):
    gtex_data = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/links/gtex_data/{tissue_name}.v8.signif_variant_gene_pairs.txt.gz', sep = '\t')
    gtex_data['chrom'] = [x.split('_')[0] for x in gtex_data['variant_id']]
    gtex_data['start'] = [int(x.split('_')[1]) for x in gtex_data['variant_id']]
    gtex_data['end'] = [int(x.split('_')[1]) for x in gtex_data['variant_id']]
    gtex_data['ref'] = [x.split('_')[2] for x in gtex_data['variant_id']]
    gtex_data['alt'] = [x.split('_')[3] for x in gtex_data['variant_id']]
    gtex_data['gene_id'] = [x.split('.')[0] for x in gtex_data['gene_id']]
    gtex_data = gtex_data[['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'maf', 'slope']]
    return(gtex_data)

def get_tss_positions(links):
    tss_starts, tss_ends = [], []
    for gene in links[['chrom_gene', 'start_gene', 'end_gene', 'strand_gene']].values:
        if gene[-1] == '+':
            tss_start = gene[1]-2500
            tss_end = gene[1]+2500
        else:
            tss_start = gene[2]-2500
            tss_end = gene[2]+2500
        tss_starts.append(tss_start)
        tss_ends.append(tss_end)
    links['start_tss_gene'] = tss_starts
    links['end_tss_gene'] = tss_ends
    return(links)

def overlap_links_and_tads(links, tads):
    links = get_tss_positions(links)
    
    overlap_links_tads_peaks = bioframe.overlap(links, tads, 
                                                cols1 = ['chrom_peak', 'start_peak', 'end_peak'], 
                                                cols2 = ['chrom_tad', 'start_tad', 'end_tad'],
                                                suffixes = ['', '_peak'], return_overlap = True)
    overlap_links_tads_genes = bioframe.overlap(overlap_links_tads_peaks, tads, 
                                                cols1 = ['chrom_gene', 'start_gene', 'end_gene'], 
                                                cols2 = ['chrom_tad', 'start_tad', 'end_tad'],
                                                suffixes = ['', '_gene'], return_overlap = True)
    overlap_links_tads = bioframe.overlap(overlap_links_tads_genes, tads, 
                                          cols1 = ['chrom_gene', 'start_tss_gene', 'end_tss_gene'], 
                                          cols2 = ['chrom_tad', 'start_tad', 'end_tad'],
                                          suffixes = ['', '_tss'], return_overlap = True)

    overlap_links_tads['tad_distance'] = np.abs(overlap_links_tads['id_tad_peak']-overlap_links_tads['id_tad_tss'])

    columns_to_drop = ['start_tss_gene', 'end_tss_gene', 
                       'chrom_tad_peak', 'start_tad_peak', 'end_tad_peak', 'id_tad_peak', 'overlap_start_peak', 'overlap_end_peak', 
                       'chrom_tad_gene', 'start_tad_gene', 'end_tad_gene', 'id_tad_gene', 'overlap_start_gene', 'overlap_end_gene',
                       'chrom_tad_tss', 'start_tad_tss', 'end_tad_tss', 'id_tad_tss', 'overlap_start_tss_gene', 'overlap_end_tss_gene']
    overlap_links_tads = overlap_links_tads.drop(columns_to_drop, axis = 1).sort_values(['link_id', 'tad_distance'])
    overlap_links_tads = overlap_links_tads.drop_duplicates(subset = overlap_links_tads.columns[:-1])
    overlap_links_tads.index = np.arange(len(overlap_links_tads))

    return(overlap_links_tads)

gtex_data_all = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/links/gtex_data/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Cells_EBV-transformed_lymphocytes.allpairs.txt.gz', sep = '\t')
gtex_data_all['chrom'] = [x.split('_')[0] for x in gtex_data_all['variant_id']]
gtex_data_all['start'] = [int(x.split('_')[1]) for x in gtex_data_all['variant_id']]
gtex_data_all['end'] = [int(x.split('_')[1]) for x in gtex_data_all['variant_id']]
gtex_data_all['ref'] = [x.split('_')[2] for x in gtex_data_all['variant_id']]
gtex_data_all['alt'] = [x.split('_')[3] for x in gtex_data_all['variant_id']]
gtex_data_all['gene_id'] = [x.split('.')[0] for x in gtex_data_all['gene_id']]
gtex_data_all = gtex_data_all[['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'maf', 'slope']]
gtex_all_genes = np.unique(gtex_data_all['gene_id'].values)

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

for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:
    print(sample)

    phased_variants = read_variants(sample)

    links_allele_specific = pd.read_csv(f'links_{sample}_allele_specific.txt', sep = '\t', header = 0)
    links_balanced = pd.read_csv(f'links_{sample}_balanced.txt', sep = '\t', header = 0)
    
    tads = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/HiC/insulation/{sample}.tads.bed', sep = '\t', header = None)
    tads.columns = ['chrom_tad', 'start_tad', 'end_tad']
    tads['id_tad'] = np.arange(len(tads))
    links_allele_specific = overlap_links_and_tads(links_allele_specific, tads)
    links_allele_specific = links_allele_specific[links_allele_specific['tad_distance'] == 0]
    links_balanced = overlap_links_and_tads(links_balanced, tads)
    links_balanced = links_balanced[links_balanced['tad_distance'] == 0]

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

    
    links_allele_specific['gene_id'] = [x.split('.')[0] for x in links_allele_specific['gene_id'].values]
    no_gene_as = len(links_allele_specific)-len(links_allele_specific[links_allele_specific['gene_id'].isin(gtex_all_genes)])
    
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
    counts_total_as = np.array([gtex_sign_count, gtex_all_count, gtex_variants_count, other_count-no_gene_as, no_gene_as])

    
    links_balanced['gene_id'] = [x.split('.')[0] for x in links_balanced['gene_id'].values]
    no_gene_balanced = len(links_balanced)-len(links_balanced[links_balanced['gene_id'].isin(gtex_all_genes)])
    
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
    counts_total_balanced = np.array([gtex_sign_count, gtex_all_count, gtex_variants_count, other_count-no_gene_balanced, no_gene_balanced])
    
    print(counts_total_as)
    print(counts_total_balanced)
    

gtex_tissues = os.listdir('/DATA/users/m.magnitov/hap_phen/links/gtex_data/')
gtex_tissues = sorted([x.split('.')[0] for x in gtex_tissues if 'all_associations' not in x and 'GTEx_Analysis' not in x])

counts_gtex_tissues = pd.DataFrame()
for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:
    print(sample)
    phased_variants = read_variants(sample)

    links_allele_specific = pd.read_csv(f'links_{sample}_allele_specific.txt', sep = '\t', header = 0)
    links_balanced = pd.read_csv(f'links_{sample}_balanced.txt', sep = '\t', header = 0)

    tads = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/HiC/insulation/{sample}.tads.bed', sep = '\t', header = None)
    tads.columns = ['chrom_tad', 'start_tad', 'end_tad']
    tads['id_tad'] = np.arange(len(tads))
    links_allele_specific = overlap_links_and_tads(links_allele_specific, tads)
    links_allele_specific = links_allele_specific[links_allele_specific['tad_distance'] == 0]
    links_balanced = overlap_links_and_tads(links_balanced, tads)
    links_balanced = links_balanced[links_balanced['tad_distance'] == 0]
    
    variants_links_allele_specific = get_variants_for_links(phased_variants, links_allele_specific)
    variants_links_balanced = get_variants_for_links(phased_variants, links_balanced)
    
    counts_gtex_tissues_as, counts_gtex_tissues_balanced = [], []
    for tissue in gtex_tissues:
        gtex_data = read_significant_gtex_data(tissue)
        variants_links_allele_specific_in_gtex = variants_links_allele_specific.merge(gtex_data, 
                                                                                      left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant', 'gene_id'], 
                                                                                      right_on = ['chrom', 'start', 'ref', 'alt', 'gene_id'])
        variants_links_balanced_in_gtex = variants_links_balanced.merge(gtex_data, 
                                                                        left_on = ['chrom_variant', 'start_variant', 'ref_variant', 'alt_variant', 'gene_id'], 
                                                                        right_on = ['chrom', 'start', 'ref', 'alt', 'gene_id'])

        counts_gtex_tissues_as.append(len(variants_links_allele_specific_in_gtex)/len(variants_links_allele_specific))
        counts_gtex_tissues_balanced.append(len(variants_links_balanced_in_gtex)/len(variants_links_balanced))
    counts_gtex_tissues[sample + '_as'] = counts_gtex_tissues_as
    counts_gtex_tissues[sample + '_balanced'] = counts_gtex_tissues_balanced
    
counts_gtex_tissues.index = gtex_tissues

eqtls_per_tissue = []
for tissue in gtex_tissues:
    gtex_data = read_significant_gtex_data(tissue)
    eqtls_per_tissue.append(len(gtex_data))
    
counts_gtex_tissues['N_eQTL'] = eqtls_per_tissue

counts_gtex_tissues.to_csv('overlap_links_same_tad_with_gtex_per_tissue.tsv', sep = '\t', header = 1, index = 1)
