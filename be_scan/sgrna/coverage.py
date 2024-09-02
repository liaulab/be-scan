"""
Author: Calvin XiaoYang Hu
Date: 240314

{Description: 
            }
"""

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import re

from _genomic_ import *
from _guideRNA_ import *

def coverage_plots(
    annotated_guides, edit_from, edit_to, protein_filepath, 

    output_name="guides_by_aa.csv", output_dir='',
    return_df=True, save_df=True,
): 
    path = Path.cwd()
    # READ GUIDES FILE #
    guides_filepath = Path(annotated_guides)
    df = pd.read_csv(guides_filepath)

    # REORGANIZE DATAFRAME BY GUIDE INTO DATAFRAME BY AMINO ACID #
    amino_acid_seq = protein_to_AAseq(protein_filepath)
    df_aa = pd.DataFrame(data=amino_acid_seq.items(), columns=['position', 'amino_acid'])

    # A LIST OF AMINO ACIDS, WHICH GUIDE IDS, GUIDE SEQUENCES, DIRECTION, MUTATIONS, MUTATION TYPES FOR EACH EDITOR TYPE #
    sgRNAID_results, seq_results, strand_results = [], [], []
    counts = []
    mut_results, muttype_results = [], []
    pre = f'{edit_from}to{edit_to}'
    for mut_pos in df_aa['amino_acid']+df_aa['position'].astype(str): 
        mut_pos = mut_pos+'[A-Z.]'
        matching_rows = df[df[f'{pre}_mutations'].str.contains(mut_pos, regex=True, na=False)]

        # ADD FINDINGS TO LISTS AND APPEND THEM TO DF #
        sgRNAIDs = matching_rows['sgRNA_ID'].tolist()
        sgRNAID_results.append(';'.join(sgRNAIDs))
        counts.append(len(sgRNAIDs))
        seqs = matching_rows['sgRNA_seq'].tolist()
        seq_results.append(';'.join(seqs))
        strands = matching_rows['sgRNA_strand'].tolist()
        strand_results.append(';'.join(strands))

        # ADD SPECIFIC SUB ENTRIES #
        muts = [x for xs in matching_rows[f'{pre}_mutations'].tolist() for x in 
                [xi for xss in xs.split(';') for xi in xss.split('/')]]
        # each list has entries that are ; and / separated so I need to split them all apart
        muts = [mut for mut in muts if bool(re.search(mut_pos, mut))]
        mut_results.append(';'.join(list(set(muts))))

    df_aa['guide_IDs'] = sgRNAID_results
    df_aa['guide_counts'] = counts
    df_aa['sgRNA_seqs'] = seq_results
    df_aa['sgRNA_strands'] = strand_results
    df_aa[f'{pre}_mutations'] = mut_results

    if save_df: 
        Path.mkdir(path / output_dir, exist_ok=True)
        df_aa.to_csv(path / output_dir / f'{pre}_{output_name}', index=False)

    # SHOW COVERAGE ALONG LENGTH OF PROTEIN PLOT AND STATS HISTOGRAM #
    plt.figure(figsize=(len(df_aa)/50, max(counts)))
    ax1 = sns.barplot(x=df_aa['position'], y=counts, 
                      color='steelblue', edgecolor='steelblue')
    ax1.set_ylabel(f"Count of {pre} Guides")
    ax1.set_title(f"Count of {pre} Guides Per Residue")
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_aa), 50), rotation = 90)

    if save_df: 
        Path.mkdir(path / output_dir, exist_ok=True)
        plt.savefig(path / output_dir / f'{pre}_{output_name.split('.')[0]}_barplot.png', dpi=400)

    # COUNT HOW MANY MUTATIONS ARE ACCESSIBLE PER RESIDUE #
    aa_list = ['R', 'H', 'K', # POS CHARGED
               'D', 'E', # NEG CHARGED
               'S', 'N', 'T', 'Q', # UNCHARGED
               'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', # HYDROPHOBIC
               'C', 'G', 'P', # SPECIAL
               '.' # TERMINATION
               ] # Y AXIS
    df_heatmap = pd.DataFrame(columns=aa_list)
    for i, (muts, aai) in enumerate(zip(mut_results, amino_acid_seq.values())): 
        muts_dict = {aa:0 for aa in aa_list}
        muts_dict[aai] = -1
        # IF THERES NO MUTATIONS #
        if len(muts) == 0: 
            df_heatmap.loc[-1] = muts_dict.values()
            continue
        # IF THERE ARE MUTATIONS #
        for mut in muts.split(';'): 
            muts_dict[mut[-1]] = 1+(2*muts_dict[mut[-1]]) # -1 STAYS -1, 0 BECOMES 1
        df_heatmap.loc[i] = muts_dict.values()
    ### its hard to distinguish between silent mutations, and when the original mutation is just there

    # SHOW WHICH MUTATIONS ARE ACCESSIBLE WITH A HEATMAP #
    plt.figure(figsize=(len(df_aa)/20, 3))
    ax2 = sns.heatmap(df_heatmap.T, cmap="vlag")
    ax2.set_xticks(range(df_heatmap.shape[0]))
    ax2.set_yticks(range(df_heatmap.shape[1]))
    ax2.set_xticklabels(df_heatmap.index, rotation=90)
    ax2.set_yticklabels(df_heatmap.columns)
    ax2.tick_params(axis='both', which='major', labelsize=4.5)
    ax2.tick_params(axis='both', which='minor', labelsize=4.5)

    if save_df: 
        Path.mkdir(path / output_dir, exist_ok=True)
        plt.savefig(path / output_dir / f'{pre}_{output_name.split('.')[0]}_heatmap.png', dpi=400)

    if return_df: return df_aa, (ax1, ax2)
