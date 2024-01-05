"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Functions for adding dataframes of guides together}
"""

from pathlib import Path
import pandas as pd

### blanks need to be filled in with blanks or No_C No_A etc
def merge_guide_df(guide_df1_filepath, guide_df2_filepath, 
                   output_name='merged_guides.csv', output_dir='', 

                   shared_col_names = ['sgRNA_seq', 'starting_frame', 'gene_pos', 
                                       'chr_pos', 'exon', 'coding_seq', 'sgRNA_strand', 
                                       'gene_strand', 'gene', 'genome_occurrences'], 
                   col1_name = 'CBE',
                   col2_name = 'ABE',

                   sort_by=['gene_pos'],
                   return_df=True, save_df=True,
                   ): 
    """
    Merge 2 Dataframes (for example ABE and CBE dataframes)
    Keeps annotated information for both separate
    Completes merge on the guide sequence, and keeps annotated information separate.
    Deletes duplicated guides between the two dataframes. 
    """
    guide_df1_filepath = Path(guide_df1_filepath)
    guide_df2_filepath = Path(guide_df2_filepath)
    df1 = pd.read_csv(guide_df1_filepath)
    df2 = pd.read_csv(guide_df2_filepath)
    assert 'sgRNA_seq' in df1.columns and 'sgRNA_seq' in df2.columns
    assert 'coding_seq' in df1.columns and 'coding_seq' in df2.columns

    shared_col_names = [name for name in df1.columns if name in df2.columns]
    # merge dataframes
    new_df = df1.merge(df2, on=shared_col_names, how='outer')
    for cond in sort_by: 
        assert cond in new_df.columns, f'Sort condition {cond} not found in both dataframe'
        new_df.sort_values(cond)

    # delete duplicates
    dupl_rows = new_df.duplicated(subset='sgRNA_seq', keep=False)
    new_df = new_df[~dupl_rows]
    dupl_rows = new_df.duplicated(subset='coding_seq', keep=False)
    new_df = new_df[~dupl_rows]

    if save_df: 
        out_filepath = Path(output_dir)
        new_df.to_csv(out_filepath / output_name, index=False)
    if return_df:
        return new_df

def add_guide_df(guides_df_filepath, additional_df_filepath,
                 output_name='new_guides.csv', output_dir='',
                 return_df=True, save_df=True,
                 ): 
    """
    Add 2 Dataframes (ex guides and control guides)
    Reassigns the columns of second dataframe into first. 
    """
    guides_df_filepath = Path(guides_df_filepath)
    additional_df_filepath = Path(additional_df_filepath)
    guides_df = pd.read_csv(guides_df_filepath)
    additional_df = pd.read_csv(additional_df_filepath)
    assert 'sgRNA_seq' in guides_df.columns and 'sgRNA_seq' in additional_df.columns

    # check column names are compatible
    assert all(name in guides_df.columns for name in additional_df.columns), "Make sure dataframe column names match."
    if 'coding_seq' not in additional_df.columns: 
        additional_df['coding_seq'] = additional_df['sgRNA_seq'] ### assuming all control guides are sense

    new_df = pd.concat([guides_df, additional_df])
    # delete duplicates
    dupl_rows = new_df.duplicated(subset='sgRNA_seq', keep=False)
    new_df = new_df[~dupl_rows]
    dupl_rows = new_df.duplicated(subset='coding_seq', keep=False)
    new_df = new_df[~dupl_rows]
    new_df = new_df[~new_df['sgRNA_seq'].isin(new_df['coding_seq']) & ~new_df['coding_seq'].isin(new_df['sgRNA_seq'])]

    if save_df: 

        new_df.to_csv(output_dir+output_name)
    if return_df:
        return new_df
