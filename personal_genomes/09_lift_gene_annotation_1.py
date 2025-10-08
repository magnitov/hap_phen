import os
import pandas as pd
import numpy as np
from gtfparse import read_gtf
from liftover import get_lifter, ChainFile

def convert_gene_annotation_row_to_haplotype(gene_annotation_instance, converting_function):
    # Convert reference coordinates to haplotype coordinates
    # In case the liftover is not successful, assign chromosome to 'error' and position to 0.

    converted_start = converting_function[gene_annotation_instance[0]][gene_annotation_instance[3]]
    converted_end = converting_function[gene_annotation_instance[0]][gene_annotation_instance[4]]
    
    if len(converted_start) != 0:
        start_chrom_hap, start_pos_hap = converted_start[0][0], converted_start[0][1]
    else:
        start_chrom_hap = 'error_chrom'
        start_pos_hap = 0
        
    if len(converted_end) != 0:
        end_chrom_hap, end_pos_hap = converted_end[0][0], converted_end[0][1]
    else:
        end_chrom_hap = 'error_chrom'
        end_pos_hap = 0

    return(start_chrom_hap, start_pos_hap, end_chrom_hap, end_pos_hap)


gene_annotation = read_gtf('/DATA/users/m.magnitov/genomes/gencode.v42.basic.annotation.gtf')
genomes_path = '/DATA/users/m.magnitov/hap_phen/personal_genomes/'

for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:

    converter_ref_to_hap1 = ChainFile(f'{genomes_path}/{sample}/genome/{sample}_ref_to_hap1.chain', 'ref', 'hap1')
    converter_ref_to_hap2 = ChainFile(f'{genomes_path}/{sample}/genome/{sample}_ref_to_hap2.chain', 'ref', 'hap2')

    lifted_gene_annotation_hap1, lifted_gene_annotation_hap2 = [], []

    for row in gene_annotation.values:

        start_chrom_hap1, start_pos_hap1, end_chrom_hap1, end_pos_hap1 = convert_gene_annotation_row_to_haplotype(row, converter_ref_to_hap1)
        lifted_gene_annotation_hap1.append([start_chrom_hap1, start_pos_hap1, end_chrom_hap1, end_pos_hap1])

        start_chrom_hap2, start_pos_hap2, end_chrom_hap2, end_pos_hap2 = convert_gene_annotation_row_to_haplotype(row, converter_ref_to_hap2)
        lifted_gene_annotation_hap2.append([start_chrom_hap2, start_pos_hap2, end_chrom_hap2, end_pos_hap2])

    lifted_gene_annotation_hap1 = pd.DataFrame(lifted_gene_annotation_hap1)
    lifted_gene_annotation_hap2 = pd.DataFrame(lifted_gene_annotation_hap2)
    lifted_gene_annotation_hap1.to_csv(f'{genomes_path}/{sample}/genome/lifted_gencode_annotation_coordinates.hap1.gtf', sep = '\t', index = 0, header = 0)
    lifted_gene_annotation_hap2.to_csv(f'{genomes_path}/{sample}/genome/lifted_gencode_annotation_coordinates.hap2.gtf', sep = '\t', index = 0, header = 0)
