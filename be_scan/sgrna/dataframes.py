

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
    
    df1 = pd.read_csv(guide_df1_filepath)
    df2 = pd.read_csv(guide_df2_filepath)

    # make sure shared_col_names appear in both dataframes
    remove_ind = []
    df1_colnames = df1.columns
    df2_colnames = df2.columns
    for i, x in enumerate(shared_col_names): 
        if x not in df1_colnames or x not in df2_colnames: 
            remove_ind.append(i)
    for i in reversed(remove_ind): 
        removed = shared_col_names.pop(i)
        print(removed, 'removed from shared_col_names')

    # gather names of other col_names and change their names according to col1_name col2_name
    df1_original_names = [x for x in df1_colnames if x not in shared_col_names]
    df1_new_names = ['_'.join([col1_name, x]) for x in df1_colnames if x not in shared_col_names]
    df2_original_names = [x for x in df2_colnames if x not in shared_col_names]
    df2_new_names = ['_'.join([col2_name, x]) for x in df2_colnames if x not in shared_col_names]
    dict1 = dict(zip(df1_original_names, df1_new_names))
    dict2 = dict(zip(df2_original_names, df2_new_names))
    df1.rename(columns=dict1, inplace=True)
    df2.rename(columns=dict2, inplace=True)

    # merge dataframes
    new_df = df1.merge(df2, on=shared_col_names, how='outer')
    for cond in sort_by: 
        assert cond in new_df.columns, f'Sort condition {cond} not found in both dataframe'
        new_df.sort_values(cond)

    # delete duplicates
    # Identify rows where either A or B is duplicated
    duplicated_rows = new_df[new_df.duplicated(subset=['sgRNA_seq'], keep=False) | new_df.duplicated(subset=['coding_seq'], keep=False)]
    # Identify rows where A appears in B
    rows_to_remove = duplicated_rows[duplicated_rows['sgRNA_seq'].isin(duplicated_rows['coding_seq'])]
    # Drop the identified rows from the original DataFrame
    new_df = new_df.drop(rows_to_remove.index)

    if save_df: 
        new_df.to_csv(output_dir+output_name)
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
    
    guides_df = pd.read_csv(guides_df_filepath)
    additional_df = pd.read_csv(additional_df_filepath)
    assert all(name in guides_df.columns for name in additional_df.columns), "Make sure dataframe column names match."

    new_df = pd.concat([guides_df, additional_df])
    
    # delete duplicates
    # Identify rows where either A or B is duplicated
    duplicated_rows = new_df[new_df.duplicated(subset=['sgRNA_seq'], keep=False) | new_df.duplicated(subset=['coding_seq'], keep=False)]
    # Identify rows where A appears in B
    rows_to_remove = duplicated_rows[duplicated_rows['sgRNA_seq'].isin(duplicated_rows['coding_seq'])]
    # Drop the identified rows from the original DataFrame
    new_df = new_df.drop(rows_to_remove.index)

    if save_df: 
        new_df.to_csv(output_dir+output_name)
    if return_df:
        return new_df
