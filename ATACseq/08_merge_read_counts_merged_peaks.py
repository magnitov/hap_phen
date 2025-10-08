import pandas as pd
import numpy as np
import functools as ft
import warnings
warnings.filterwarnings('ignore')

def get_counts_table(peaks, sample, replicate, counts_path):
    counts_hap1 = pd.read_csv(f'{counts_path}/counts_all_peaks_{sample}_{replicate}_hap1.txt', 
                              sep = '\t', header = None, names = ['chrom', 'start', 'end', 'peak_id', f'{sample}_{replicate}_hap1'])
    counts_hap2 = pd.read_csv(f'{counts_path}/counts_all_peaks_{sample}_{replicate}_hap2.txt', 
                              sep = '\t', header = None, names = ['chrom', 'start', 'end', 'peak_id', f'{sample}_{replicate}_hap2'])
    
    peaks = peaks.merge(counts_hap1, on = 'peak_id', how = 'left')
    peaks = peaks.merge(counts_hap2, on = 'peak_id', how = 'left')
    
    peaks = peaks.fillna(0)
    
    peaks = peaks[['seqname', 'start', 'end', 'peak_id', f'{sample}_{replicate}_hap1', f'{sample}_{replicate}_hap2']]
    
    return(peaks)


counts_path = '/DATA/users/m.magnitov/hap_phen/ATACseq/counts'
peaks_path = '/DATA/users/m.magnitov/hap_phen/ATACseq/peaks'

for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:
    # Read peaks
    peaks = pd.read_csv(peaks_path + '/all_peaks.canonical.replicated.no_blacklist.bed',
                        sep = '\t', header = None, names = ['seqname', 'start', 'end', 'peak_id'])

    # Read peak counts
    counts_allelic_rep1 = get_counts_table(peaks, sample, 'rep1', counts_path)
    counts_allelic_rep2 = get_counts_table(peaks, sample, 'rep2', counts_path)
    counts_allelic_rep3 = get_counts_table(peaks, sample, 'rep3', counts_path)
    counts_unassigned_rep1 = get_counts_table(peaks, sample, 'rep1_unassigned', counts_path)
    counts_unassigned_rep2 = get_counts_table(peaks, sample, 'rep2_unassigned', counts_path)
    counts_unassigned_rep3 = get_counts_table(peaks, sample, 'rep3_unassigned', counts_path)

    # Create a table with allelic counts
    count_tables = [counts_allelic_rep1, counts_allelic_rep2, counts_allelic_rep3]
    counts_allelic = ft.reduce(lambda left, right: pd.merge(left, right, on = ['seqname', 'start', 'end', 'peak_id']), count_tables)
    for column in counts_allelic.columns[4:]:
        counts_allelic[column] = [int(x) for x in counts_allelic[column]]
    counts_allelic.drop(counts_allelic.columns[[0, 1, 2]], axis = 1).to_csv(f'{counts_path}/counts_all_peaks_allelic_{sample}.txt', sep = '\t', header = 1, index = 0)

    # Create a table with unassigned counts
    count_tables = [counts_unassigned_rep1, counts_unassigned_rep2, counts_unassigned_rep3]
    counts_unassigned = ft.reduce(lambda left, right: pd.merge(left, right, on = ['seqname', 'start', 'end', 'peak_id']), count_tables)
    for column in counts_unassigned.columns[4:]:
        counts_unassigned[column] = [int(x) for x in counts_unassigned[column]]
    counts_unassigned.drop(counts_unassigned.columns[[0, 1, 2]], axis = 1).to_csv(f'{counts_path}/counts_all_peaks_unassigned_{sample}.txt', sep = '\t', header = 1, index = 0)

    # Create a table with total counts per peak
    counts_total = pd.DataFrame()
    counts_total['peak_id'] = counts_allelic['peak_id'].values
    counts_total[sample + '_rep1'] = (counts_unassigned[sample + '_rep1_unassigned_hap1'] + counts_unassigned[sample + '_rep1_unassigned_hap2'])/2 + \
                                      counts_allelic[sample + '_rep1_hap1'] + counts_allelic[sample + '_rep1_hap2']
    counts_total[sample + '_rep2'] = (counts_unassigned[sample + '_rep2_unassigned_hap1'] + counts_unassigned[sample + '_rep2_unassigned_hap2'])/2 + \
                                      counts_allelic[sample + '_rep2_hap1'] + counts_allelic[sample + '_rep2_hap2']
    counts_total[sample + '_rep3'] = (counts_unassigned[sample + '_rep3_unassigned_hap1'] + counts_unassigned[sample + '_rep3_unassigned_hap2'])/2 + \
                                      counts_allelic[sample + '_rep3_hap1'] + counts_allelic[sample + '_rep3_hap2']
    counts_total[sample + '_rep1'] = [int(x) for x in counts_total[sample + '_rep1']]
    counts_total[sample + '_rep2'] = [int(x) for x in counts_total[sample + '_rep2']]
    counts_total[sample + '_rep3'] = [int(x) for x in counts_total[sample + '_rep3']]
    counts_total.to_csv(f'{counts_path}/counts_all_peaks_total_{sample}.txt', sep = '\t', header = 1, index = 0)

