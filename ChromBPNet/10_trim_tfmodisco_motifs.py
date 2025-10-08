import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgridspec
import h5py as h5
from genome_tools.plotting import sequence
import warnings
warnings.filterwarnings('ignore')

def extract_pwm(sample, pwm_id):
    with open(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/tf_modisco/{sample}_hg38.modisco_results.meme') as f:
        pwms = f.read()

    pwm_to_extract = pwms.split(pwm_id)[1].split('MOTIF')[0].split('\n')
    pwm_to_extract = [x.lstrip(' ').split() for x in pwm_to_extract if x != ''][1:]
    pwm_to_extract = pd.DataFrame(pwm_to_extract)
    for i in range(4):
        pwm_to_extract[i] = [float(x) for x in pwm_to_extract[i]]
    pwm_to_extract = pwm_to_extract.values
    return(pwm_to_extract)

def relative_info_content(pwm):
    p = pwm/np.sum(pwm, axis = 1)[:,np.newaxis]
    ic = 2+np.sum(p*np.nan_to_num(np.log2(p)), axis = 1)
    ric = p*ic[:,np.newaxis]
    return(ric)

def plot_pwm(pwm):
    w = pwm.shape[0]

    fig = plt.figure()
    fig.set_size_inches(w*0.125, 0.5)

    figw, figh = fig.get_size_inches() 

    gs = mgridspec.GridSpec(1, 1)
    gs.update(left=0, right=1, top=1, bottom=0)

    ax = fig.add_subplot(gs[:,:])

    sequence.seq_plot(relative_info_content(pwm), ax=ax)

    ax.set_xlim(left=0, right=w)
    ax.set_ylim(bottom=0, top=2.1)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    [ax.spines[loc].set_visible(False) for loc in ['top', 'bottom', 'left', 'right']]
    
for sample in ['NA12878', 'NA18983', 'HG01241', 'HG02601', 'HG03464']:
    os.mkdir(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/tf_modisco/motifs_{sample}_hg38')
    fh = h5.File(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/tf_modisco/{sample}_hg38.modisco_results.h5')
    patterns = list(fh['pos_patterns'])
    for pattern in patterns:
        # read pwm
        pwm = extract_pwm(sample, pattern)
        #plot_pwm(pwm)

        # trim pwm
        ic = relative_info_content(pwm)
        total_ic = ic.sum(axis = 1)
        cdf = np.cumsum(total_ic)/np.sum(total_ic)
        s = np.where(cdf > 0.05)[0][0]
        e = np.where(cdf > 0.95)[0][0] + 1    
        pwm_trimmed = pwm[s:e,:]
        #plot_pwm(pwm_trimmed)
        
        # save trimmed pwm
        header_line = f'{sample}_{pattern}'
        mat = pd.DataFrame(pwm_trimmed.T, index=['A:', 'C:', 'G:', 'T:']).to_string(header = False)
        with open(f'/DATA/users/m.magnitov/hap_phen/chromBPNet/tf_modisco/motifs_{sample}_hg38/{sample}_{pattern}.pfm', 'w') as fh:
            fh.write(header_line + '\n' + mat + '\n')
