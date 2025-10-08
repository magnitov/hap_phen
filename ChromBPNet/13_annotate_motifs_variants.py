import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from Bio import SeqIO
import bioframe
from multiprocessing import Pool
import traceback
import warnings
warnings.filterwarnings('ignore')


###
### Functions to load and prepare data
###
def read_variants(sample):
    variants = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/{sample}.phased.het.bed', 
                           sep = '\t', header = None)
    variants[11] = [int(x.split('|')[0]) for x in variants[10].values]
    variants[12] = [int(x.split('|')[1].split(':')[0]) for x in variants[10].values]
    variants = variants[[0, 2, 2, 3, 5, 6, 11, 12]]
    variants.columns = ['chrom', 'start', 'end', 'variant_id', 'ref', 'alt', 'allele_ref', 'allele_alt']
    return(variants)


def read_variants_on_haplotype(sample, haplotype):
    with open(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/genome/{sample}_{haplotype}.vcf', 'r') as v_hap:
        variants_hap = v_hap.readlines()
    variants_hap = [x.rstrip('\n').split('\t') for x in variants_hap if x[0] != '#']
    variants_hap = pd.DataFrame(variants_hap)
    variants_hap = variants_hap[[0, 1, 1, 2]]
    variants_hap.columns = ['chrom', 'start_' + haplotype, 'end_' + haplotype, 'variant_id']
    return(variants_hap)


def read_asocr_peaks(sample):
    asocr_peaks = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/ATACseq/asocr/{sample}_allele_specific.bed', sep = '\t', header = None)
    asocr_peaks.columns = ['chrom', 'start', 'end', 'peak_id', 'baseMean', 'log2FoldChange', 'FDR']
    return(asocr_peaks[['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR']])


def read_balanced_peaks(sample):
    balanced_peaks = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/ATACseq/asocr/{sample}_balanced.bed', sep = '\t', header = None)
    balanced_peaks.columns = ['chrom', 'start', 'end', 'peak_id', 'baseMean', 'log2FoldChange', 'FDR']
    return(balanced_peaks[['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR']])


def get_variants_for_peaks(phased_variants, peaks):
    peaks = peaks[peaks['chrom'] != 'chrX']
    phased_variants_overlap_peaks = bioframe.overlap(phased_variants, peaks, suffixes = ('', '_peak'))
    phased_variants_overlap_peaks = phased_variants_overlap_peaks.dropna()

    putative_non_overlaps = peaks[~peaks['peak_id'].isin(phased_variants_overlap_peaks['peak_id_peak'])]['peak_id'].values
    if len(putative_non_overlaps) != 0:
        peaks_non_overlaps = peaks[peaks['peak_id'].isin(putative_non_overlaps)]

        phased_variants_overlap_non_overlaps = bioframe.closest(peaks_non_overlaps, phased_variants, suffixes = ('_peak', ''))
        sort_order = list(phased_variants_overlap_non_overlaps.columns[5:]) + list(phased_variants_overlap_non_overlaps.columns[:5])
        phased_variants_overlap_non_overlaps = phased_variants_overlap_non_overlaps[sort_order]
        phased_variants_overlap_non_overlaps = phased_variants_overlap_non_overlaps.drop('distance', axis = 1)

    phased_variants_overlap_peaks = pd.concat([phased_variants_overlap_peaks, phased_variants_overlap_non_overlaps], ignore_index = True)
    phased_variants_overlap_peaks = phased_variants_overlap_peaks.sort_values(['chrom', 'start'])
    phased_variants_overlap_peaks.index = np.arange(len(phased_variants_overlap_peaks))
    return(phased_variants_overlap_peaks)

###
### Functions to extract data about each haplotype
###
def get_variant_coordinates_in_haplotypes(phased_variants_data, variant_id):
    chrom_name = phased_variants_data[phased_variants_data['variant_id'] == variant_id]['chrom'].values[0]
    pos_hap1 = int(phased_variants_data[phased_variants_data['variant_id'] == variant_id][f'start_hap1'].values[0])
    pos_hap2 = int(phased_variants_data[phased_variants_data['variant_id'] == variant_id][f'start_hap2'].values[0])
    return(chrom_name, pos_hap1, pos_hap2)

def extract_sequence_for_region(genome, chrom, pos, flank, variant_margin):
    seq = str(genome[chrom][pos-flank-1:pos+flank+variant_margin])
    return(seq)

def extract_contribution_scores_for_region(contrib_scores, chrom, pos, flank, variant_margin):
    contrib_scores_per_region = contrib_scores[(contrib_scores['chrom'] == chrom) & (contrib_scores['start'] >= pos-flank-1) & (contrib_scores['end'] <= pos+flank+variant_margin)]
    #print(contrib_scores_per_region, len(contrib_scores_per_region))
    if len(contrib_scores_per_region) == 0:
        # in some cases contribution scores are just zeros, when variant is outside of peak summit
        contrib_zeros = pd.DataFrame({'chrom': [chrom]*(2*flank+variant_margin+1),
                                      'start': np.arange(pos-flank-1, pos+flank+variant_margin),
                                      'end': np.arange(pos-flank, pos+flank+variant_margin+1)})
        contrib_zeros.columns = ['chrom', 'start', 'end']
        contrib_zeros['start'] = [int(x) for x in contrib_zeros['start']]
        contrib_zeros['end'] = [int(x) for x in contrib_zeros['end']]
        if 'score_hap1' in contrib_scores_per_region.columns:
            contrib_zeros['score_hap1'] = np.zeros(len(contrib_zeros))
        else:
            contrib_zeros['score_hap2'] = np.zeros(len(contrib_zeros))
        contrib_zeros['relative_pos'] = np.arange(pos-flank, pos+flank+variant_margin+1)-pos
        return(contrib_zeros)
    
    elif len(contrib_scores_per_region) != 2*flank+variant_margin+1:
        # in some cases contribution scores are not binned per 1 bp, but are merged across 2 bp
        # this has to be split into 1 bp windows
        contrib_corrected = pd.DataFrame({'chrom': [chrom]*(2*flank+variant_margin+1),
                                          'start': np.arange(pos-flank-1, pos+flank+variant_margin),
                                          'end': np.arange(pos-flank, pos+flank+variant_margin+1)})
        contrib_corrected.columns = ['chrom', 'start', 'end']
        contrib_corrected['start'] = [int(x) for x in contrib_corrected['start']]
        contrib_corrected['end'] = [int(x) for x in contrib_corrected['end']]

        contrib_corrected = bioframe.overlap(contrib_corrected, contrib_scores_per_region, how='left', suffixes=('_new','_old'))
        #print(contrib_corrected, len(contrib_corrected))
        if 'score_hap1' in contrib_scores_per_region.columns:
            contrib_corrected = contrib_corrected[['chrom_new', 'start_new', 'end_new', 'score_hap1_old']]
            contrib_corrected.columns = ['chrom', 'start', 'end', 'score_hap1']
        else:
            contrib_corrected = contrib_corrected[['chrom_new', 'start_new', 'end_new', 'score_hap2_old']]
            contrib_corrected.columns = ['chrom', 'start', 'end', 'score_hap2']
        #print(contrib_corrected, len(contrib_corrected))
        contrib_corrected['relative_pos'] = np.arange(pos-flank, pos+flank+variant_margin+1)-pos
        return(contrib_corrected)
    
    else:
        contrib_scores_per_region['relative_pos'] = contrib_scores_per_region['end']-pos
        return(contrib_scores_per_region)
    
def contribution_data_to_matrix(variant_contrib_data):
    # Create matrices by taking into account Ns for InDels
    contrib_score_matrix_with_n_hap1, contrib_score_matrix_with_n_hap2 = [], []
    for i in range(len(variant_contrib_data)):
        if variant_contrib_data['seq_hap1'].values[i] == 'A':
            contrib_score_matrix_with_n_hap1.append([variant_contrib_data['score_hap1'].values[i], 0, 0, 0])
        elif variant_contrib_data['seq_hap1'].values[i] == 'C':
            contrib_score_matrix_with_n_hap1.append([0, variant_contrib_data['score_hap1'].values[i], 0, 0])
        elif variant_contrib_data['seq_hap1'].values[i] == 'G':
            contrib_score_matrix_with_n_hap1.append([0, 0, variant_contrib_data['score_hap1'].values[i], 0])
        elif variant_contrib_data['seq_hap1'].values[i] == 'T':
            contrib_score_matrix_with_n_hap1.append([0, 0, 0, variant_contrib_data['score_hap1'].values[i]])
        else:
            contrib_score_matrix_with_n_hap1.append([0, 0, 0, 0])

    for i in range(len(variant_contrib_data)):
        if variant_contrib_data['seq_hap2'].values[i] == 'A':
            contrib_score_matrix_with_n_hap2.append([variant_contrib_data['score_hap2'].values[i], 0, 0, 0])
        elif variant_contrib_data['seq_hap2'].values[i] == 'C':
            contrib_score_matrix_with_n_hap2.append([0, variant_contrib_data['score_hap2'].values[i], 0, 0])
        elif variant_contrib_data['seq_hap2'].values[i] == 'G':
            contrib_score_matrix_with_n_hap2.append([0, 0, variant_contrib_data['score_hap2'].values[i], 0])
        elif variant_contrib_data['seq_hap2'].values[i] == 'T':
            contrib_score_matrix_with_n_hap2.append([0, 0, 0, variant_contrib_data['score_hap2'].values[i]])
        else:
            contrib_score_matrix_with_n_hap2.append([0, 0, 0, 0])
            
    # Create matrices by not taking into account Ns, matrices for hap1 and hap2 will be of different lengths in this case
    contrib_score_matrix_hap1, contrib_score_matrix_hap2 = [], []
    for i in range(len(variant_contrib_data)):
        if variant_contrib_data['seq_hap1'].values[i] == 'A':
            contrib_score_matrix_hap1.append([variant_contrib_data['score_hap1'].values[i], 0, 0, 0])
        elif variant_contrib_data['seq_hap1'].values[i] == 'C':
            contrib_score_matrix_hap1.append([0, variant_contrib_data['score_hap1'].values[i], 0, 0])
        elif variant_contrib_data['seq_hap1'].values[i] == 'G':
            contrib_score_matrix_hap1.append([0, 0, variant_contrib_data['score_hap1'].values[i], 0])
        elif variant_contrib_data['seq_hap1'].values[i] == 'T':
            contrib_score_matrix_hap1.append([0, 0, 0, variant_contrib_data['score_hap1'].values[i]])
        else:
            pass

    for i in range(len(variant_contrib_data)):
        if variant_contrib_data['seq_hap2'].values[i] == 'A':
            contrib_score_matrix_hap2.append([variant_contrib_data['score_hap2'].values[i], 0, 0, 0])
        elif variant_contrib_data['seq_hap2'].values[i] == 'C':
            contrib_score_matrix_hap2.append([0, variant_contrib_data['score_hap2'].values[i], 0, 0])
        elif variant_contrib_data['seq_hap2'].values[i] == 'G':
            contrib_score_matrix_hap2.append([0, 0, variant_contrib_data['score_hap2'].values[i], 0])
        elif variant_contrib_data['seq_hap2'].values[i] == 'T':
            contrib_score_matrix_hap2.append([0, 0, 0, variant_contrib_data['score_hap2'].values[i]])
        else:
            pass

    contrib_score_matrix_with_n_hap1 = pd.DataFrame(contrib_score_matrix_with_n_hap1)
    contrib_score_matrix_with_n_hap2 = pd.DataFrame(contrib_score_matrix_with_n_hap2)
    contrib_score_matrix_with_n_hap1.columns = ['A', 'C', 'G', 'T']
    contrib_score_matrix_with_n_hap2.columns = ['A', 'C', 'G', 'T']
    contrib_score_matrix_with_n_hap1['relative_pos'] = variant_contrib_data['relative_pos']
    contrib_score_matrix_with_n_hap2['relative_pos'] = variant_contrib_data['relative_pos']
    
    contrib_score_matrix_hap1 = pd.DataFrame(contrib_score_matrix_hap1)
    contrib_score_matrix_hap2 = pd.DataFrame(contrib_score_matrix_hap2)
    contrib_score_matrix_hap1.columns = ['A', 'C', 'G', 'T']
    contrib_score_matrix_hap2.columns = ['A', 'C', 'G', 'T']
    contrib_score_matrix_hap1['relative_pos'] = variant_contrib_data['relative_pos']
    contrib_score_matrix_hap2['relative_pos'] = variant_contrib_data['relative_pos']

    return(contrib_score_matrix_with_n_hap1, contrib_score_matrix_with_n_hap2, contrib_score_matrix_hap1, contrib_score_matrix_hap2)

def read_pwm(pwm_name):
    pwm = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/tf_modisco_pwms/{pwm_name}.pfm', 
                      header = None, skiprows = 1, sep = '\s+')
    pwm.index = [x.rstrip(':') for x in pwm[0].values]
    pwm = pwm.drop([0], axis = 1).T
    return(pwm)

def annotate_variant_with_motifs(variant_id, phased_variants_data, hap1_genome, hap2_genome, contrib_scores_hap1, contrib_scores_hap2, 
                                 flank = 25, correlation_threshold = 0.7):
    print(variant_id)
    
    # get variant margin in case of indel
    variant_phasing_data = phased_variants_data[phased_variants_data['variant_id'] == variant_id][['ref', 'alt', 'allele_ref', 'allele_alt']]
    variant_phasing_data = list(variant_phasing_data.values[0])

    if variant_phasing_data[2] == 0:
        variant_allele_hap1 = variant_phasing_data[0]
        variant_allele_hap2 = variant_phasing_data[1]
    else:
        variant_allele_hap1 = variant_phasing_data[1]
        variant_allele_hap2 = variant_phasing_data[0]

    variant_margin_hap1 = len(variant_allele_hap1)-1
    variant_margin_hap2 = len(variant_allele_hap2)-1
    
    if variant_margin_hap1 == 0 and variant_margin_hap2 == 0:
        variant_type = 'snp'
    else:
        variant_type = 'indel'

    # get sequences and contribution scores of both haplotypes
    chrom_name, variant_pos_hap1, variant_pos_hap2 = get_variant_coordinates_in_haplotypes(phased_variants_data, variant_id)
    seq_hap1 = extract_sequence_for_region(hap1_genome, chrom_name, variant_pos_hap1, flank, variant_margin_hap1)
    seq_hap2 = extract_sequence_for_region(hap2_genome, chrom_name, variant_pos_hap2, flank, variant_margin_hap2)
    contrib_hap1 = extract_contribution_scores_for_region(contrib_scores_hap1, chrom_name, variant_pos_hap1, flank, variant_margin_hap1)
    contrib_hap2 = extract_contribution_scores_for_region(contrib_scores_hap2, chrom_name, variant_pos_hap2, flank, variant_margin_hap2)

    # align contribution scores in case of indel
    if len(contrib_hap1)>len(contrib_hap2):
        updated_relative_positions = np.hstack((contrib_hap2['relative_pos'].values[:flank+1], contrib_hap2['relative_pos'].values[flank+1:]+variant_margin_hap1))
        contrib_hap2['relative_pos'] = updated_relative_positions
        seq_hap2 = seq_hap2[:flank+1] + 'N'*variant_margin_hap1 + seq_hap2[flank+1:]
    elif len(contrib_hap1)<len(contrib_hap2):
        updated_relative_positions = np.hstack((contrib_hap1['relative_pos'].values[:flank+1], contrib_hap1['relative_pos'].values[flank+1:]+variant_margin_hap2))
        contrib_hap1['relative_pos'] = updated_relative_positions
        seq_hap1 = seq_hap1[:flank+1] + 'N'*variant_margin_hap2 + seq_hap1[flank+1:]
    else:
        pass

    # merge sequences and contribution scores into one dataframe
    contrib = contrib_hap1.merge(contrib_hap2, on = ['chrom', 'relative_pos'], how = 'outer').sort_values('relative_pos').fillna(0)
    seq = pd.DataFrame({'relative_pos': np.arange(min(contrib['relative_pos']), max(contrib['relative_pos'])+1), 
                        'seq_hap1': list(seq_hap1), 'seq_hap2': list(seq_hap2)})
    variant_contrib_data = contrib.merge(seq, on = 'relative_pos')
    variant_contrib_data = variant_contrib_data[['relative_pos', 'score_hap1', 'score_hap2', 'seq_hap1', 'seq_hap2']]
    variant_contrib_matrix_hap1, variant_contrib_matrix_hap2, variant_contrib_matrix_without_n_hap1, variant_contrib_matrix_without_n_hap2 = contribution_data_to_matrix(variant_contrib_data)

    # calculate correlation of motifs with contribution scores per haplotype
    correlations_patterns = []
    pwm_name_list = os.listdir('/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/tf_modisco_pwms/')
    pwm_name_list = sorted([x.split('.')[0].replace("/", "-") for x in pwm_name_list if 'pfm' in x])

    variant_margin = max(variant_margin_hap1, variant_margin_hap2)
    for pattern in pwm_name_list:
        
        # extract pattern matrix
        pattern_matrix = read_pwm(pattern)
        
        # trim contribution score matrices based on pattern size and variant margin
        variant_contrib_matrix_trimmed_hap1 = variant_contrib_matrix_hap1[(variant_contrib_matrix_hap1['relative_pos'] > -len(pattern_matrix)) &\
                                                                          (variant_contrib_matrix_hap1['relative_pos'] < len(pattern_matrix) + variant_margin)]
        variant_contrib_matrix_trimmed_hap2 = variant_contrib_matrix_hap2[(variant_contrib_matrix_hap2['relative_pos'] > -len(pattern_matrix)) &\
                                                                          (variant_contrib_matrix_hap2['relative_pos'] < len(pattern_matrix) + variant_margin)]

        # sweep the pattern matrix across the contribution score matrices
        for i in range(len(pattern_matrix) + variant_margin):

            pattern_hap1_to_compare = variant_contrib_matrix_trimmed_hap1[i:i+len(pattern_matrix)].drop(['relative_pos'], axis = 1)
            pattern_hap2_to_compare = variant_contrib_matrix_trimmed_hap2[i:i+len(pattern_matrix)].drop(['relative_pos'], axis = 1)
            
            # evaluate in regular direction
            corr_hap1 = pearsonr(pattern_hap1_to_compare.values.flatten(), pattern_matrix.values.flatten())[0]
            corr_hap2 = pearsonr(pattern_hap2_to_compare.values.flatten(), pattern_matrix.values.flatten())[0]
            
            # evaluate in reverse complement direction
            corr_hap1_rc = pearsonr(pattern_hap1_to_compare.values.flatten(), pattern_matrix.values[::-1,::-1].flatten())[0]
            corr_hap2_rc = pearsonr(pattern_hap2_to_compare.values.flatten(), pattern_matrix.values[::-1,::-1].flatten())[0]
            
            # combine into dataframe
            correlations_patterns.append([variant_contrib_matrix_trimmed_hap1[i:i+len(pattern_matrix)]['relative_pos'].values[0], 
                                          corr_hap1, corr_hap2, pattern, len(pattern_matrix), 0,
                                          np.sum(np.abs(pattern_hap1_to_compare.values.flatten())), 
                                          np.sum(np.abs(pattern_hap2_to_compare.values.flatten()))])
            correlations_patterns.append([variant_contrib_matrix_trimmed_hap1[i:i+len(pattern_matrix)]['relative_pos'].values[0], 
                                          corr_hap1_rc, corr_hap2_rc, pattern, len(pattern_matrix), 1,
                                          np.sum(np.abs(pattern_hap1_to_compare.values.flatten())), 
                                          np.sum(np.abs(pattern_hap2_to_compare.values.flatten()))])
            
        if variant_type == 'indel':
            # trim contribution score matrices based on pattern size and variant margin
            variant_contrib_matrix_trimmed_hap1 = variant_contrib_matrix_without_n_hap1[(variant_contrib_matrix_without_n_hap1['relative_pos'] > -len(pattern_matrix)) &\
                                                                                        (variant_contrib_matrix_without_n_hap1['relative_pos'] < len(pattern_matrix))]
            variant_contrib_matrix_trimmed_hap2 = variant_contrib_matrix_without_n_hap2[(variant_contrib_matrix_without_n_hap2['relative_pos'] > -len(pattern_matrix)) &\
                                                                                        (variant_contrib_matrix_without_n_hap2['relative_pos'] < len(pattern_matrix))]
            # sweep the pattern matrix across the contribution score matrices
            for i in range(len(pattern_matrix)):

                pattern_hap1_to_compare = variant_contrib_matrix_trimmed_hap1[i:i+len(pattern_matrix)].drop(['relative_pos'], axis = 1)
                pattern_hap2_to_compare = variant_contrib_matrix_trimmed_hap2[i:i+len(pattern_matrix)].drop(['relative_pos'], axis = 1)

                # evaluate in regular direction
                corr_hap1 = pearsonr(pattern_hap1_to_compare.values.flatten(), pattern_matrix.values.flatten())[0]
                corr_hap2 = pearsonr(pattern_hap2_to_compare.values.flatten(), pattern_matrix.values.flatten())[0]

                # evaluate in reverse complement direction
                corr_hap1_rc = pearsonr(pattern_hap1_to_compare.values.flatten(), pattern_matrix.values[::-1,::-1].flatten())[0]
                corr_hap2_rc = pearsonr(pattern_hap2_to_compare.values.flatten(), pattern_matrix.values[::-1,::-1].flatten())[0]

                # combine into dataframe
                correlations_patterns.append([variant_contrib_matrix_trimmed_hap1[i:i+len(pattern_matrix)]['relative_pos'].values[0], 
                                              corr_hap1, corr_hap2, pattern, len(pattern_matrix), 0,
                                              np.sum(np.abs(pattern_hap1_to_compare.values.flatten())), 
                                              np.sum(np.abs(pattern_hap2_to_compare.values.flatten()))])
                correlations_patterns.append([variant_contrib_matrix_trimmed_hap1[i:i+len(pattern_matrix)]['relative_pos'].values[0], 
                                              corr_hap1_rc, corr_hap2_rc, pattern, len(pattern_matrix), 1,
                                              np.sum(np.abs(pattern_hap1_to_compare.values.flatten())), 
                                              np.sum(np.abs(pattern_hap2_to_compare.values.flatten()))])

    # create dataframe with correlation scores
    correlations_patterns = pd.DataFrame(correlations_patterns).fillna(0)
    correlations_patterns.columns = ['relative_pos_start', 'r_hap1', 'r_hap2', 'pattern_id', 'pattern_length', 'is_rc', 'cs_motif_hap1', 'cs_motif_hap2']
    
    # filter dataframe by concordance of correlation and contribution score difference
    correlations_patterns['delta_cs_motif'] = correlations_patterns['cs_motif_hap2'] - correlations_patterns['cs_motif_hap1']
    correlations_patterns['delta_r'] = correlations_patterns['r_hap2'] - correlations_patterns['r_hap1']
    correlations_patterns = correlations_patterns[(correlations_patterns['delta_cs_motif'] * correlations_patterns['delta_r']) > 0]
    
    # filter dataframe by correlation and return the results
    correlations_patterns['r_max'] = [np.nanmax([x, y]) for (x, y) in correlations_patterns[['r_hap1', 'r_hap2']].values]
    correlations_patterns = correlations_patterns.sort_values(['r_max'], ascending = False)
    correlations_patterns = correlations_patterns.drop_duplicates(subset = ['pattern_id', 'relative_pos_start'])
    correlations_patterns = correlations_patterns[correlations_patterns['r_max'] > correlation_threshold]
    correlations_patterns = correlations_patterns.drop_duplicates(subset = ['pattern_id', 'relative_pos_start'])
    correlations_patterns['rank'] = np.arange(len(correlations_patterns))+1
    
    # add variant data for output
    correlations_patterns['chrom'] = [chrom_name]*max(len(correlations_patterns), 1)
    correlations_patterns['start'] = [variant_id.split(':')[1]]*len(correlations_patterns)
    correlations_patterns['variant_id'] = [variant_id]*len(correlations_patterns)
    correlations_patterns['start_hap1'] = [variant_pos_hap1]*len(correlations_patterns)
    correlations_patterns['start_hap2'] = [variant_pos_hap2]*len(correlations_patterns)
    correlations_patterns['cs_flank_hap1'] = [sum(variant_contrib_data['score_hap1'])]*len(correlations_patterns)
    correlations_patterns['cs_flank_hap2'] = [sum(variant_contrib_data['score_hap2'])]*len(correlations_patterns)
    correlations_patterns['delta_cs_flank'] = correlations_patterns['cs_flank_hap2']-correlations_patterns['cs_flank_hap1']
    
    correlations_patterns = correlations_patterns[['chrom', 'start', 'variant_id', 
                                                   'pattern_id', 'is_rc', 'rank', 'relative_pos_start', 'pattern_length', 
                                                   'r_hap1', 'r_hap2', 'r_max', 'delta_r', 
                                                   'cs_motif_hap1', 'cs_motif_hap2', 'delta_cs_motif',
                                                   'cs_flank_hap1', 'cs_flank_hap2', 'delta_cs_flank',
                                                   'start_hap1', 'start_hap2']]
    
    return(correlations_patterns)

def args_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield(lst[i:i + n])

#############################################################################

###
### Read data for sample
###
sample = 'HG03464'
print('Reading data')

# Read phased variants coordinates in hg38 and exclude chrX from analyses
phased_variants = read_variants(sample)
phased_variants = phased_variants[phased_variants['chrom'] != 'chrX']

# Read phased variants coordinates lifted to haplotype 1 and haplotype 2 coordinates
variants_hap1 = read_variants_on_haplotype(sample, 'hap1')
variants_hap2 = read_variants_on_haplotype(sample, 'hap2')

# Read allele-specific and balanced accessibility peaks
asocr_peaks = read_asocr_peaks(sample)
balanced_peaks = read_balanced_peaks(sample)

# Overlap phased variants with allele-specific peaks
phased_variants_overlap_asocr = get_variants_for_peaks(phased_variants, asocr_peaks)
phased_variants_overlap_asocr = phased_variants_overlap_asocr.merge(variants_hap1, on = ['chrom', 'variant_id'])
phased_variants_overlap_asocr = phased_variants_overlap_asocr.merge(variants_hap2, on = ['chrom', 'variant_id'])

# Overlap phased variants with balanced peaks
phased_variants_overlap_balanced = get_variants_for_peaks(phased_variants, balanced_peaks)
phased_variants_overlap_balanced = phased_variants_overlap_balanced.merge(variants_hap1, on = ['chrom', 'variant_id'])
phased_variants_overlap_balanced = phased_variants_overlap_balanced.merge(variants_hap2, on = ['chrom', 'variant_id'])

# Read genomes for each haplotype
hap1_genome, hap2_genome = {}, {}
for record in SeqIO.parse(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/genome/{sample}_hap1.fa', 'fasta'):
    hap1_genome[record.id] = str(record.seq)
for record in SeqIO.parse(f'/DATA/users/m.magnitov/hap_phen/personal_genomes/{sample}/genome/{sample}_hap2.fa', 'fasta'):
    hap2_genome[record.id] = str(record.seq)
    
# Read contribution scores for each haplotype
contrib_scores_hap1 = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/contrib/{sample}_hap1.profile_scores.bed', sep = '\t', header = None)
contrib_scores_hap2 = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/contrib/{sample}_hap2.profile_scores.bed', sep = '\t', header = None)
contrib_scores_hap1.columns = ['chrom', 'start', 'end', 'score_hap1']
contrib_scores_hap2.columns = ['chrom', 'start', 'end', 'score_hap2']

###
### Annotate variants with motifs
###
print('Annotating data')

#
# Variants in allele-specific peaks
#
#args_asocr = [(phased_variants_overlap_asocr['variant_id'].values[i], phased_variants_overlap_asocr, 
#               hap1_genome, hap2_genome, contrib_scores_hap1, contrib_scores_hap2) for i in range(len(phased_variants_overlap_asocr))]        
#args_asocr = list(args_chunks(args_asocr, 200))
#
#for batch in range(0, len(args_asocr)):
#    if __name__ == '__main__':
#        p = Pool(processes = 8)
#        variants_annotation_motifs_asocr = p.starmap(annotate_variant_with_motifs, args_asocr[batch])
#
#    print('Saving annotated ALLELE-SPECIFIC variants')
#    variants_annotation_motifs_asocr = pd.DataFrame(np.concatenate(variants_annotation_motifs_asocr))
#    variants_annotation_motifs_asocr.columns = ['chrom', 'start', 'variant_id', 
#                                                'pattern_id', 'is_rc', 'rank', 'relative_pos_start', 'pattern_length', 
#                                                'r_hap1', 'r_hap2', 'r_max', 'delta_r', 
#                                                'cs_hap1', 'cs_hap2', 'delta_cs',
#                                                'cs_flank_hap1', 'cs_flank_hap2', 'delta_cs_flank',
#                                                'start_hap1', 'start_hap2']
#    variants_annotation_motifs_asocr.to_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/{sample}_variants_annotation_asocr_batch{batch}.txt', sep = '\t', header = 1, index = 0)
               
#
# Variants in balanced peaks
#
args_balanced = [(phased_variants_overlap_balanced['variant_id'].values[i], phased_variants_overlap_balanced, 
               hap1_genome, hap2_genome, contrib_scores_hap1, contrib_scores_hap2) for i in range(len(phased_variants_overlap_balanced))]        
args_balanced = list(args_chunks(args_balanced, 200))

for batch in range(0, len(args_balanced)):
    if __name__ == '__main__':
        p = Pool(processes = 8)
        variants_annotation_motifs_balanced = p.starmap(annotate_variant_with_motifs, args_balanced[batch])

    print('Saving annotated BALANCED variants')
    variants_annotation_motifs_balanced = pd.DataFrame(np.concatenate(variants_annotation_motifs_balanced))
    variants_annotation_motifs_balanced.columns = ['chrom', 'start', 'variant_id', 
                                                'pattern_id', 'is_rc', 'rank', 'relative_pos_start', 'pattern_length', 
                                                'r_hap1', 'r_hap2', 'r_max', 'delta_r', 
                                                'cs_hap1', 'cs_hap2', 'delta_cs',
                                                'cs_flank_hap1', 'cs_flank_hap2', 'delta_cs_flank',
                                                'start_hap1', 'start_hap2']
    variants_annotation_motifs_balanced.to_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/{sample}_variants_annotation_balanced_batch{batch}.txt', sep = '\t', header = 1, index = 0)
