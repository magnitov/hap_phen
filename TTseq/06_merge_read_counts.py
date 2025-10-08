import pandas as pd
import numpy as np
from gtfparse import read_gtf
import functools as ft
import warnings
warnings.filterwarnings('ignore')

def get_counts_table(genes, sample, replicate, counts_path):
    counts_hap1_fwd = pd.read_csv(f'{counts_path}/counts_{sample}_{replicate}_hap1_forward.txt', 
                                  sep = '\t', header = None, names = ['gene_id', f'{sample}_{replicate}_hap1_fwd'])
    counts_hap1_rev = pd.read_csv(f'{counts_path}/counts_{sample}_{replicate}_hap1_reverse.txt', 
                                  sep = '\t', header = None, names = ['gene_id', f'{sample}_{replicate}_hap1_rev'])
    counts_hap2_fwd = pd.read_csv(f'{counts_path}/counts_{sample}_{replicate}_hap2_forward.txt', 
                                  sep = '\t', header = None, names = ['gene_id', f'{sample}_{replicate}_hap2_fwd'])
    counts_hap2_rev = pd.read_csv(f'{counts_path}/counts_{sample}_{replicate}_hap2_reverse.txt', 
                                  sep = '\t', header = None, names = ['gene_id', f'{sample}_{replicate}_hap2_rev'])
    
    genes = genes.merge(counts_hap1_fwd, on = 'gene_id', how = 'left')
    genes = genes.merge(counts_hap1_rev, on = 'gene_id', how = 'left')
    genes = genes.merge(counts_hap2_fwd, on = 'gene_id', how = 'left')
    genes = genes.merge(counts_hap2_rev, on = 'gene_id', how = 'left')
    
    counts_hap1, counts_hap2 = [], []
    for gene in genes.values:
        if gene[3] == '+':
            counts_hap1.append(gene[7])
            counts_hap2.append(gene[9])
        else:
            counts_hap1.append(gene[8])
            counts_hap2.append(gene[10])
    genes[f'{sample}_{replicate}_hap1'] = counts_hap1
    genes[f'{sample}_{replicate}_hap2'] = counts_hap2
    
    genes = genes.fillna(0)
    
    genes = genes[['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name', f'{sample}_{replicate}_hap1', f'{sample}_{replicate}_hap2']]
    
    return(genes)

# Read GTF with gene annotation
genes = read_gtf('/DATA/users/m.magnitov/genomes/gencode.v42.basic.annotation.gtf')
genes = genes[genes['feature'] == 'gene']
genes = genes[['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name']]

counts_path = '/DATA/users/m.magnitov/hap_phen/TTseq/counts'

for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:
    # Read gene counts
    counts_allelic_rep1 = get_counts_table(genes, sample, 'rep1', counts_path)
    counts_allelic_rep2 = get_counts_table(genes, sample, 'rep2', counts_path)
    counts_unassigned_rep1 = get_counts_table(genes, sample, 'rep1_unassigned', counts_path)
    counts_unassigned_rep2 = get_counts_table(genes, sample, 'rep2_unassigned', counts_path)

    # Create a table with allelic counts
    count_tables = [counts_allelic_rep1, counts_allelic_rep2]
    counts_allelic = ft.reduce(lambda left, right: pd.merge(left, right, on = ['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name']), count_tables)
    for column in counts_allelic.columns[7:]:
        counts_allelic[column] = [int(x) for x in counts_allelic[column]]
    counts_allelic.drop(counts_allelic.columns[[0, 1, 2, 3, 5, 6]], axis = 1).to_csv(f'{counts_path}/counts_allelic_{sample}.txt', sep = '\t', header = 1, index = 0)

    # Create a table with unassigned counts
    count_tables = [counts_unassigned_rep1, counts_unassigned_rep2]
    counts_unassigned = ft.reduce(lambda left, right: pd.merge(left, right, on = ['seqname', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name']), count_tables)
    for column in counts_unassigned.columns[7:]:
        counts_unassigned[column] = [int(x) for x in counts_unassigned[column]]
    counts_unassigned.drop(counts_unassigned.columns[[0, 1, 2, 3, 5, 6]], axis = 1).to_csv(f'{counts_path}/counts_unassigned_{sample}.txt', sep = '\t', header = 1, index = 0)

    # Create a table with total counts per gene
    counts_total = pd.DataFrame()
    counts_total['gene_id'] = counts_allelic['gene_id'].values
    counts_total[sample + '_rep1'] = (counts_unassigned[sample + '_rep1_unassigned_hap1'] + counts_unassigned[sample + '_rep1_unassigned_hap2'])/2 + \
                                      counts_allelic[sample + '_rep1_hap1'] + counts_allelic[sample + '_rep1_hap2']
    counts_total[sample + '_rep2'] = (counts_unassigned[sample + '_rep2_unassigned_hap1'] + counts_unassigned[sample + '_rep2_unassigned_hap2'])/2 + \
                                      counts_allelic[sample + '_rep2_hap1'] + counts_allelic[sample + '_rep2_hap2']
    counts_total[sample + '_rep1'] = [int(x) for x in counts_total[sample + '_rep1']]
    counts_total[sample + '_rep2'] = [int(x) for x in counts_total[sample + '_rep2']]
    counts_total.to_csv(f'{counts_path}/counts_total_{sample}.txt', sep = '\t', header = 1, index = 0)

