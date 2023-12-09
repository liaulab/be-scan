"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Run count_reads, merge_and_norm, average_reps, compare_conds}
"""

from be_scan.analysis.count_reads import count_reads
from be_scan.analysis.merge_and_norm import merge_and_norm
from be_scan.analysis.average_reps import average_reps
from be_scan.analysis.compare_conds import compare_conds

def batch_process(sample_sheet, in_ref, in_comparisons, 
                  
                  file_dir='', save=True, return_df=False,

                  KEY_INTERVAL=(10,80), KEY='CGAAACACC', KEY_REV='GTTTTAGA', dont_trim_G=False,

                  t0='t0', dir_counts='', 
                  out_reads='agg_log2_t0.csv', 
                  out_conds='agg_t0_conds.csv', 

                  out_comps='agg_comps.csv', 
                 ):
    
    """
    Given a set of sgRNA sequences, a sample sheet, and a list of comparisons, 
    complete count_reads, merge_and_norm, average_reps, compare_conds runs.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'fastq_file', 'counts_file', 'noncounts_file', 'stats_file'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_np (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    in_comparisons : in_comparisons .csv in format (name, treatment, control)
        A dataframe denoting the comparisons to make, with the comparison
        being treatment - control. The output column
        headers will be labeled by the name in the dataframe.

    file_dir : str or path, defaults to ''
        String or path to the directory where all files are found and saved. 
    save : bool, default True
        Choose files for export to csv. The default is 'all', which is the
        aggregated read counts ('reads'), log2 normalized values ('log2'), and
        t0 normalized values ('t0'). You may also enter any combination of
        'reads', 'log2', 't0') as a list of strings to choose which ones to
        save ('all' is equivalent to a list of all three). None will not export
        any files to csv.
    return_df : {None, 'reads', 'log2', 't0'}, default None
        Whether to return a dataframe at function end. The default is None,
        which returns nothing. However, you can return the reads, log2 norm,
        or t0 norm values by calling 'reads', 'log2', or't0', respectively.

    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_REV : str, default 'GTTTTAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.

    t0 : str, default 't0'
        Name of the t0 sample in dict_counts. If you have multiple t0 samples
        (e.g. paired t0s for specific samples), then you will need to run this
        function separately for each set of samples with their appropriate t0.
    dir_counts : str, default ''
        Name of the subfolder to find the read count csv files. The default is
        the current working directory.
    out_reads : str, default 'agg_reads.csv'
        Name of the aggregated raw reads csv output file.
    out_log2 : str, default 'agg_log2.csv'
        Name of the aggregated log2 normalized values csv output file.
    out_t0 : str, default 'agg_t0_reps.csv'
        Name of the aggregated t0 normalized values csv output file.

    out_conds : str, default 'agg_t0_conds.csv'
        Name of the averaged replicate values csv output file.

    out_comps : str, default 'agg_comps.csv'
        Name of the comparisons csv output file.
    """

    count_reads(sample_sheet=sample_sheet, 
                in_ref=in_ref, 
                file_dir=file_dir,

                KEY_INTERVAL=KEY_INTERVAL, 
                KEY=KEY, 
                KEY_REV=KEY_REV, 
                dont_trim_G=dont_trim_G,
                )
    
    merge_and_norm(sample_sheet=sample_sheet, 
                   in_ref=file_dir+'counts_library.csv', 
                   file_dir=file_dir,

                   save=save, 
                   out=out_reads, 
                   return_df=return_df,
                   )
    
    average_reps(sample_sheet=sample_sheet, 
                 in_lfc=file_dir+out_reads, 
                 file_dir=file_dir,

                 save=save, 
                 out_conds=out_conds, 
                 return_df=return_df,
                 )
    
    compare_conds(in_comparisons=in_comparisons, 
                  in_conds=file_dir+out_conds, 
                  file_dir=file_dir,

                  out_comps=out_comps, 
                  save=save,
                  return_df=return_df,
                  )
    
