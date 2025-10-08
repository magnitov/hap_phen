import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from Bio import SeqIO
import bioframe
from multiprocessing import Pool
import warnings
warnings.filterwarnings('ignore')


###
### Functions to load and prepare data
###
def read_asocr_peaks(sample):
    asocr_peaks = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/ATACseq/asocr/{sample}_allele_specific.bed', sep = '\t', header = None)
    asocr_peaks.columns = ['chrom', 'start', 'end', 'peak_id', 'baseMean', 'log2FoldChange', 'FDR']
    return(asocr_peaks[['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR']])

def read_balanced_peaks(sample):
    balanced_peaks = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/ATACseq/asocr/{sample}_balanced.bed', sep = '\t', header = None)
    balanced_peaks.columns = ['chrom', 'start', 'end', 'peak_id', 'baseMean', 'log2FoldChange', 'FDR']
    return(balanced_peaks[['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR']])

###
### Functions to extract data about each haplotype
###
def extract_sequence_for_region(genome, chrom, start, end):
    return(str(genome[chrom][start:end]))

def extract_contribution_scores_for_region(contrib_scores, chrom, start, end):
    contrib_scores_per_region = contrib_scores[(contrib_scores['chrom'] == chrom) & (contrib_scores['start'] >= start) & (contrib_scores['end'] <= end)]
    #print(contrib_scores_per_region, len(contrib_scores_per_region))
    if len(contrib_scores_per_region) == 0:
        # in some cases contribution scores are just zeros, when variant is outside of peak summit
        contrib_zeros = pd.DataFrame({'chrom': [chrom]*(end-start+1),
                                      'start': np.arange(start, end),
                                      'end': np.arange(start+1, end+1)})
        contrib_zeros.columns = ['chrom', 'start', 'end']
        contrib_zeros['start'] = [int(x) for x in contrib_zeros['start']]
        contrib_zeros['end'] = [int(x) for x in contrib_zeros['end']]
        if 'score_hap1' in contrib_scores_per_region.columns:
            contrib_zeros['score_hap1'] = np.zeros(len(contrib_zeros))
        else:
            contrib_zeros['score_hap2'] = np.zeros(len(contrib_zeros))
        contrib_zeros['relative_pos'] = np.arange(len(contrib_zeros))
        return(contrib_zeros)
    
    elif len(contrib_scores_per_region) != end-start:
        # in some cases contribution scores are not binned per 1 bp, but are merged across 2 bp
        # this has to be split into 1 bp windows
        contrib_corrected = pd.DataFrame({'chrom': [chrom]*(end-start),
                                          'start': np.arange(start, end),
                                          'end': np.arange(start+1, end+1)})
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
        contrib_corrected['relative_pos'] = np.arange(len(contrib_corrected))
        return(contrib_corrected)
    
    else:
        contrib_scores_per_region['relative_pos'] = np.arange(len(contrib_scores_per_region))
        return(contrib_scores_per_region)

def contribution_data_to_matrix(variant_contrib_data):
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
            contrib_score_matrix_hap1.append([0, 0, 0, 0])

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
            contrib_score_matrix_hap2.append([0, 0, 0, 0])

    contrib_score_matrix_hap1 = pd.DataFrame(contrib_score_matrix_hap1)
    contrib_score_matrix_hap2 = pd.DataFrame(contrib_score_matrix_hap2)
    contrib_score_matrix_hap1.columns = ['A', 'C', 'G', 'T']
    contrib_score_matrix_hap2.columns = ['A', 'C', 'G', 'T']
    contrib_score_matrix_hap1['relative_pos'] = variant_contrib_data['relative_pos']
    contrib_score_matrix_hap2['relative_pos'] = variant_contrib_data['relative_pos']

    return(contrib_score_matrix_hap1, contrib_score_matrix_hap2)

def read_pwm(pwm_name):
    pwm = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/tf_modisco_pwms/{pwm_name}.pfm', 
                      header = None, skiprows = 1, sep = '\s+')
    pwm.index = [x.rstrip(':') for x in pwm[0].values]
    pwm = pwm.drop([0], axis = 1).T
    return(pwm)

def annotate_peak_with_motifs(peak_id, peaks_data, contrib_scores_hap1, contrib_scores_hap2, correlation_threshold = 0.7, contribution_threshold = 0.01):
    print(peak_id)
    
    # get sequences and contribution scores of both haplotypes
    chrom_peak, start_peak, end_peak = peaks_data[peaks_data['peak_id'] == peak_id].values[0][:3]
    seq_hg38 = extract_sequence_for_region(hg38_genome, chrom_peak, start_peak, end_peak)
    contrib_hap1 = extract_contribution_scores_for_region(contrib_scores_hap1, chrom_peak, start_peak, end_peak)
    contrib_hap2 = extract_contribution_scores_for_region(contrib_scores_hap2, chrom_peak, start_peak, end_peak)

    # merge sequences and contribution scores into one dataframe
    contrib = contrib_hap1.merge(contrib_hap2, on = ['chrom', 'relative_pos'], how = 'outer').sort_values('relative_pos').fillna(0)
    seq = pd.DataFrame({'relative_pos': np.arange(len(seq_hg38)), 'seq_hap1': list(seq_hg38), 'seq_hap2': list(seq_hg38)})
    variant_contrib_data = contrib.merge(seq, on = 'relative_pos')
    variant_contrib_data = variant_contrib_data[['relative_pos', 'score_hap1', 'score_hap2', 'seq_hap1', 'seq_hap2']]
    variant_contrib_matrix_hap1, variant_contrib_matrix_hap2 = contribution_data_to_matrix(variant_contrib_data)
    
    # calculate correlation of motifs with contribution scores per haplotype
    correlations_patterns = []
    pwm_name_list = os.listdir('/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/tf_modisco_pwms/')
    pwm_name_list = sorted([x.split('.')[0].replace("/", "-") for x in pwm_name_list if 'pfm' in x])

    for pattern in pwm_name_list:
        # extract pattern matrix
        pattern_matrix = read_pwm(pattern)

        # sweep the pattern matrix across the contribution score matrices
        for i in range(len(variant_contrib_data) - len(pattern_matrix)):

            pattern_hap1_to_compare = variant_contrib_matrix_hap1[i:i+len(pattern_matrix)].drop(['relative_pos'], axis = 1)
            pattern_hap2_to_compare = variant_contrib_matrix_hap2[i:i+len(pattern_matrix)].drop(['relative_pos'], axis = 1)

            # evaluate in regular direction
            corr_hap1 = pearsonr(pattern_hap1_to_compare.values.flatten(), pattern_matrix.values.flatten())[0]
            corr_hap2 = pearsonr(pattern_hap2_to_compare.values.flatten(), pattern_matrix.values.flatten())[0]
            # evaluate in reverse complement direction
            corr_hap1_rc = pearsonr(pattern_hap1_to_compare.values.flatten(), pattern_matrix.values[::-1,::-1].flatten())[0]
            corr_hap2_rc = pearsonr(pattern_hap2_to_compare.values.flatten(), pattern_matrix.values[::-1,::-1].flatten())[0]

            # combine into dataframe
            correlations_patterns.append([variant_contrib_matrix_hap1[i:i+len(pattern_matrix)]['relative_pos'].values[0], 
                                          corr_hap1, corr_hap2, pattern, len(pattern_matrix), 0,
                                          np.sum(np.abs(pattern_hap1_to_compare.values.flatten())), 
                                          np.sum(np.abs(pattern_hap2_to_compare.values.flatten()))])
            correlations_patterns.append([variant_contrib_matrix_hap2[i:i+len(pattern_matrix)]['relative_pos'].values[0], 
                                          corr_hap1_rc, corr_hap2_rc, pattern, len(pattern_matrix), 1,
                                          np.sum(np.abs(pattern_hap1_to_compare.values.flatten())), 
                                          np.sum(np.abs(pattern_hap2_to_compare.values.flatten()))])

    # create dataframe with correlation scores
    correlations_patterns = pd.DataFrame(correlations_patterns)
    correlations_patterns.columns = ['relative_pos_start', 'r_hap1', 'r_hap2', 'pattern_id', 'pattern_length', 'is_rc', 'cs_hap1', 'cs_hap2']
    
    # filter dataframe by concordance of correlation and contribution score difference
    correlations_patterns['delta_cs'] = correlations_patterns['cs_hap2'] - correlations_patterns['cs_hap1']
    correlations_patterns['delta_r'] = correlations_patterns['r_hap2'] - correlations_patterns['r_hap1']
    correlations_patterns = correlations_patterns[(correlations_patterns['delta_cs'] * correlations_patterns['delta_r']) > 0]
    
    # filter dataframe by correlation and return the results
    correlations_patterns['r_max'] = [max(x, y) for (x, y) in correlations_patterns[['r_hap1', 'r_hap2']].values]
    correlations_patterns = correlations_patterns[(correlations_patterns['r_hap1'] > correlation_threshold) &\
                                                  (correlations_patterns['r_hap2'] > correlation_threshold)]
    
    # filter dataframe by constibutions and return the results
    correlations_patterns['cs_max'] = [max(x, y) for (x, y) in correlations_patterns[['cs_hap1', 'cs_hap2']].values]
    correlations_patterns = correlations_patterns[(correlations_patterns['cs_hap1'] > contribution_threshold) &\
                                                  (correlations_patterns['cs_hap2'] > contribution_threshold)]
    
    # add variant data for output
    correlations_patterns['chrom'] = [chrom_peak]*max(len(correlations_patterns), 1)
    correlations_patterns['start'] = [start_peak]*len(correlations_patterns)
    correlations_patterns['end'] = [end_peak]*len(correlations_patterns)
    correlations_patterns['peak_id'] = [peak_id]*len(correlations_patterns)
    
    correlations_patterns = correlations_patterns[['chrom', 'start', 'end', 'peak_id', 
                                                   'pattern_id', 'is_rc', 'relative_pos_start', 'pattern_length', 
                                                   'r_hap1', 'r_hap2', 'delta_r', 
                                                   'cs_hap1', 'cs_hap2', 'delta_cs']].sort_values('relative_pos_start')
    
    return(correlations_patterns)

def args_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield(lst[i:i + n])
        
#############################################################################

###
### Read data for sample
###
sample = 'NA12878'
print('Reading data')

# Read allele-specific and balanced accessibility peaks
asocr_peaks = read_asocr_peaks(sample)
balanced_peaks = read_balanced_peaks(sample)
asocr_peaks = asocr_peaks[asocr_peaks['chrom'] != 'chrX']
balanced_peaks = balanced_peaks[balanced_peaks['chrom'] != 'chrX']

peaks_hap1 = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/predict_peaks/{sample}_hap1.narrowPeak', sep = '\t', header = None)
peaks_hap2 = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/predict_peaks/{sample}_hap2.narrowPeak', sep = '\t', header = None)
peaks_all = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/ATACseq/peaks/{sample}_merged_peaks.narrowPeak', sep = '\t', header = None)

asocr_peaks = asocr_peaks.merge(peaks_all, left_on = ['chrom', 'start', 'end'], right_on = [0, 1, 2]).merge(peaks_hap1, left_on = [3], right_on = [3]).merge(peaks_hap2, left_on = [3], right_on = [3])
asocr_peaks = asocr_peaks[['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR', '1_y', '2_y', 1, 2]]
asocr_peaks.columns = ['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR', 'start_hap1', 'end_hap1', 'start_hap2', 'end_hap2']

balanced_peaks = balanced_peaks.merge(peaks_all, left_on = ['chrom', 'start', 'end'], right_on = [0, 1, 2]).merge(peaks_hap1, left_on = [3], right_on = [3]).merge(peaks_hap2, left_on = [3], right_on = [3])
balanced_peaks = balanced_peaks[['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR', '1_y', '2_y', 1, 2]]
balanced_peaks.columns = ['chrom', 'start', 'end', 'peak_id', 'log2FoldChange', 'FDR', 'start_hap1', 'end_hap1', 'start_hap2', 'end_hap2']

# Read hg38 genome
hg38_genome = {}
for record in SeqIO.parse('/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa', 'fasta'):
    hg38_genome[record.id] = str(record.seq)
    
# Read contribution scores for each haplotype
contrib_scores_hap1 = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/contrib/{sample}_hap1.profile_scores.hg38.bed', sep = '\t', header = None)
contrib_scores_hap2 = pd.read_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/contrib/{sample}_hap2.profile_scores.hg38.bed', sep = '\t', header = None)
contrib_scores_hap1.columns = ['chrom', 'start', 'end', 'score_hap1']
contrib_scores_hap2.columns = ['chrom', 'start', 'end', 'score_hap2']

###
### Motifs in allele-specific peaks
###
print('Annotating data')

args_asocr = [(asocr_peaks['peak_id'].values[i], asocr_peaks, contrib_scores_hap1, contrib_scores_hap2) for i in range(len(asocr_peaks))]
args_asocr = list(args_chunks(args_asocr, 200))

for batch in range(0, len(args_asocr)):
    if __name__ == '__main__':
        p = Pool(processes = 8)
        peaks_annotation_motifs_asocr = p.starmap(annotate_peak_with_motifs, args_asocr[batch])

    print('Saving annotated ALLELE-SPECIFIC peaks')
    peaks_annotation_motifs_asocr = pd.DataFrame(np.concatenate(peaks_annotation_motifs_asocr))
    peaks_annotation_motifs_asocr.columns = ['chrom', 'start', 'end', 'peak_id', 
                                             'pattern_id', 'is_rc', 'relative_pos_start', 'pattern_length', 
                                             'r_hap1', 'r_hap2', 'delta_r', 
                                             'cs_hap1', 'cs_hap2', 'delta_cs']
    peaks_annotation_motifs_asocr.to_csv(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/variants_annotation/{sample}_peaks_annotation_asocr_batch{batch}.txt', sep = '\t', header = 1, index = 0)
    