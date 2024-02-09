import argparse
import inspect

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main():

### analysis ###

    parser = argparse.ArgumentParser(prog="be_scan")
    subparsers = parser.add_subparsers(required=True)

    from be_scan.analysis import validate_cloning
    signat_vc = inspect.signature(validate_cloning)
    parser_validate_cloning = subparsers.add_parser('validate_cloning',
                                                    help=next(line for line in validate_cloning.__doc__.splitlines() if line),
                                                    description=validate_cloning.__doc__,
                                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_validate_cloning.add_argument('query_dir', type=str)
    parser_validate_cloning.add_argument('spacers', type=str)
    parser_validate_cloning.add_argument('vector', type=str)
    parser_validate_cloning.add_argument('out_csv', type=str, help="Output CSV file name")
    parser_validate_cloning.add_argument('--enzyme', type=str, default="Esp3I")
    parser_validate_cloning.add_argument('--overhangs', type=str, default=('CACC', 'GTTT'), nargs=2)
    parser_validate_cloning.add_argument('--ref_name_pattern', type=str, default=signat_vc.parameters['ref_name_pattern'].default)
    parser_validate_cloning.add_argument('--flank_width', type=int, default=signat_vc.parameters['flank_width'].default)
    parser_validate_cloning.set_defaults(func=validate_cloning)

    ##################################################

    from be_scan.analysis import count_reads
    signat_cr = inspect.signature(count_reads)
    parser_count_reads = subparsers.add_parser('count_reads',
                                               help=next(line for line in count_reads.__doc__.splitlines() if line),
                                               description=count_reads.__doc__,
                                               formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_count_reads.add_argument('sample_sheet', type=str)
    parser_count_reads.add_argument('annotated_lib', type=str)
    parser_count_reads.add_argument('--file_dir', type=str, default=signat_cr.parameters['file_dir'].default)
    parser_count_reads.add_argument('--KEY_INTERVAL', type=int, nargs=2, default=signat_cr.parameters['KEY_INTERVAL'].default)
    parser_count_reads.add_argument('--KEY', type=str, default=signat_cr.parameters['KEY'].default)
    parser_count_reads.add_argument('--KEY_REV', type=str, default=signat_cr.parameters['KEY_REV'].default)
    parser_count_reads.add_argument('--dont_trim_G', action='store_true')
    parser_count_reads.add_argument('--out_dir', type=str, default=signat_cr.parameters['out_dir'].default)
    parser_count_reads.add_argument('--out_file', type=str, default=signat_cr.parameters['out_file'].default)
    parser_count_reads.add_argument('--save', action='store_false', default=signat_cr.parameters['save'].default)
    parser_count_reads.add_argument('--return_df', action='store_false', default=signat_cr.parameters['return_df'].default)
    parser_count_reads.add_argument('--save_files', action='store_false', default=signat_cr.parameters['save_files'].default)
    parser_count_reads.set_defaults(func=count_reads)

    ##################################################
    
    from be_scan.analysis import merge_and_norm
    signat_man = inspect.signature(merge_and_norm)
    parser_merge_and_norm = subparsers.add_parser('merge_and_norm', 
                                                  help=next(line for line in merge_and_norm.__doc__.splitlines() if line),
                                                  description=merge_and_norm.__doc__,
                                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_merge_and_norm.add_argument('sample_sheet', type=str)
    parser_merge_and_norm.add_argument('annotated_lib', type=str)
    parser_merge_and_norm.add_argument('--controls', nargs='+', type=str, default=signat_man.parameters['controls'].default)
    parser_merge_and_norm.add_argument('--out_dir', type=str, default=signat_man.parameters['out_dir'].default)
    parser_merge_and_norm.add_argument('--out_file', type=str, default=signat_man.parameters['out_file'].default)
    parser_merge_and_norm.add_argument('--save', action='store_false', default=signat_man.parameters['save'].default)
    parser_merge_and_norm.add_argument('--return_df', action='store_false', default=signat_man.parameters['return_df'].default)
    parser_merge_and_norm.set_defaults(func=merge_and_norm)
    
    ##################################################
    
    from be_scan.analysis import average_reps
    signat_ar = inspect.signature(average_reps)
    parser_average_reps = subparsers.add_parser('average_reps', 
                                                help=next(line for line in average_reps.__doc__.splitlines() if line),
                                                description=average_reps.__doc__,
                                                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_average_reps.add_argument('sample_sheet', type=str)
    parser_average_reps.add_argument('annotated_lib', type=str)
    parser_average_reps.add_argument('--out_dir', type=str, default=signat_ar.parameters['out_dir'].default)
    parser_average_reps.add_argument('--out_file', type=str, default=signat_ar.parameters['out_file'].default)
    parser_average_reps.add_argument('--save', action='store_false', default=signat_ar.parameters['save'].default)
    parser_average_reps.add_argument('--return_df', action='store_false', default=signat_ar.parameters['return_df'].default)
    parser_average_reps.set_defaults(func=average_reps)

    ##################################################
    
    from be_scan.analysis import compare_conds
    signat_cc = inspect.signature(compare_conds)
    parser_compare_conds = subparsers.add_parser('compare_conds', 
                                                help=next(line for line in compare_conds.__doc__.splitlines() if line),
                                                description=compare_conds.__doc__,
                                                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_compare_conds.add_argument('comparisons', type=str)
    parser_compare_conds.add_argument('annotated_lib', type=str)
    parser_compare_conds.add_argument('--out_dir', type=str, default=signat_cc.parameters['out_dir'].default)
    parser_compare_conds.add_argument('--out_file', type=str, default=signat_cc.parameters['out_file'].default)
    parser_compare_conds.add_argument('--save', action='store_false', default=signat_cc.parameters['save'].default)
    parser_compare_conds.add_argument('--return_df', action='store_false', default=signat_cc.parameters['return_df'].default)
    parser_compare_conds.set_defaults(func=compare_conds)
    
    ##################################################
    
    from be_scan.analysis import calc_controls
    signat_ccs = inspect.signature(calc_controls)
    parser_calc_controls = subparsers.add_parser('calc_controls', 
                                                help=next(line for line in calc_controls.__doc__.splitlines() if line),
                                                description=calc_controls.__doc__,
                                                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_calc_controls.add_argument('annotated_lib', type=str)
    parser_calc_controls.add_argument('neg_ctrl_col', type=str)
    parser_calc_controls.add_argument('-sc', '--stats_comparisons', nargs='+', type=str, required=True)
    parser_calc_controls.add_argument('-ncc', '--neg_ctrl_conditions', nargs='+', type=str, required=True)
    parser_calc_controls.add_argument('--out_dir', type=str, default=signat_ccs.parameters['out_dir'].default)
    parser_calc_controls.add_argument('--out_file', type=str, default=signat_ccs.parameters['out_file'].default)
    parser_calc_controls.add_argument('--save', action='store_false', default=signat_ccs.parameters['save'].default)
    parser_calc_controls.add_argument('--return_txt', action='store_false', default=signat_ccs.parameters['return_txt'].default)
    parser_calc_controls.set_defaults(func=calc_controls)
    
    ##################################################
    
    from be_scan.analysis import batch_process
    signat_bp = inspect.signature(batch_process)
    parser_batch_process = subparsers.add_parser('batch_process', 
                                                 help=next(line for line in batch_process.__doc__.splitlines() if line),
                                                 description=batch_process.__doc__,
                                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_batch_process.add_argument('sample_sheet', type=str)
    parser_batch_process.add_argument('annotated_lib', type=str)
    parser_batch_process.add_argument('comparisons', type=str)
    parser_batch_process.add_argument('neg_ctrl_col', type=str)
    parser_batch_process.add_argument('-ncc', '--neg_ctrl_conditions', nargs='+', type=str, required=True)
    
    parser_batch_process.add_argument('--KEY_INTERVAL', type=int, nargs=2, default=signat_bp.parameters['KEY_INTERVAL'].default)
    parser_batch_process.add_argument('--KEY', type=str, default=signat_bp.parameters['KEY'].default)
    parser_batch_process.add_argument('--KEY_REV', type=str, default=signat_bp.parameters['KEY_REV'].default)
    parser_batch_process.add_argument('--dont_trim_G', action='store_true')
    parser_batch_process.add_argument('--file_dir', type=str, default=signat_bp.parameters['file_dir'].default)
    parser_batch_process.add_argument('--controls', type=str, default=signat_bp.parameters['controls'].default)

    parser_batch_process.add_argument('--out_dir', type=str, default=signat_bp.parameters['out_dir'].default)
    parser_batch_process.add_argument('--out_counts', type=str, default=signat_bp.parameters['out_counts'].default) # count_reads
    parser_batch_process.add_argument('--out_lfc', type=str, default=signat_bp.parameters['out_lfc'].default) # merge_and_norm
    parser_batch_process.add_argument('--out_conds', type=str, default=signat_bp.parameters['out_conds'].default) # average_reps
    parser_batch_process.add_argument('--out_comps', type=str, default=signat_bp.parameters['out_comps'].default) # compare_conds
    parser_batch_process.add_argument('--out_stats', type=str, default=signat_bp.parameters['out_stats'].default) # calc_controls
    parser_batch_process.add_argument('--save', action='store_false', default=signat_bp.parameters['save'].default)
    parser_batch_process.add_argument('--return_df', action='store_true', default=signat_bp.parameters['return_df'].default)
    parser_batch_process.set_defaults(func=batch_process)
    
    ##################################################

### plot ###

    from be_scan.plot import plot_scatterplot
    signat_ps = inspect.signature(plot_scatterplot)
    parser_plot_scatterplot = subparsers.add_parser('plot_scatterplot', 
                                                    help=next(line for line in plot_scatterplot.__doc__.splitlines() if line),
                                                    description=plot_scatterplot.__doc__,
                                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_plot_scatterplot.add_argument('df_filepath', type=str)
    parser_plot_scatterplot.add_argument('x_column', type=str)
    parser_plot_scatterplot.add_argument('-c', '--comparisons', nargs='+', type=str, required=True)
    # optional variables
    parser_plot_scatterplot.add_argument('--filter_val', action='store_true', default=signat_ps.parameters['filter_val'].default)
    parser_plot_scatterplot.add_argument('--val_cols', nargs='+', type=str, default=signat_ps.parameters['val_cols'].default)
    parser_plot_scatterplot.add_argument('--val_min', type=float, default=signat_ps.parameters['val_min'].default)
    parser_plot_scatterplot.add_argument('--filter_params', action='store_true', default=signat_ps.parameters['filter_params'].default)
    parser_plot_scatterplot.add_argument('--params_cols', nargs='+', type=str, default=signat_ps.parameters['params_cols'].default)
    parser_plot_scatterplot.add_argument('--params_conditions', nargs='+', type=list, default=signat_ps.parameters['params_conditions'].default)
    parser_plot_scatterplot.add_argument('--include_hue', action='store_true', default=signat_ps.parameters['include_hue'].default)
    parser_plot_scatterplot.add_argument('--hue_col', type=str, default=signat_ps.parameters['hue_col'].default)
    parser_plot_scatterplot.add_argument('--hue_order', nargs='+', type=str, default=signat_ps.parameters['hue_order'].default)
    parser_plot_scatterplot.add_argument('--palette', nargs='+', type=str, default=signat_ps.parameters['palette'].default)
    parser_plot_scatterplot.add_argument('--neg_ctrl', action='store_true', default=signat_ps.parameters['neg_ctrl'].default)
    parser_plot_scatterplot.add_argument('--neg_ctrl_col', type=str, default=signat_ps.parameters['neg_ctrl_col'].default)
    parser_plot_scatterplot.add_argument('--neg_ctrl_conditions', nargs='+', type=str, default=signat_ps.parameters['neg_ctrl_conditions'].default)
    parser_plot_scatterplot.add_argument('--autoannot', action='store_true', default=signat_ps.parameters['autoannot'].default)
    parser_plot_scatterplot.add_argument('--autoannot_label', type=str, default=signat_ps.parameters['autoannot_label'].default)
    parser_plot_scatterplot.add_argument('--autoannot_top', type=int, default=signat_ps.parameters['autoannot_top'].default)
    parser_plot_scatterplot.add_argument('--autoannot_cutoff', type=float, default=signat_ps.parameters['autoannot_cutoff'].default)
    parser_plot_scatterplot.add_argument('--xlab', type=str, default=signat_ps.parameters['xlab'].default)
    parser_plot_scatterplot.add_argument('--ylab', type=str, default=signat_ps.parameters['ylab'].default)
    parser_plot_scatterplot.add_argument('--col_label', type=str, default=signat_ps.parameters['col_label'].default)
    parser_plot_scatterplot.add_argument('--savefig', action='store_false', default=signat_ps.parameters['savefig'].default)
    parser_plot_scatterplot.add_argument('--out_name', type=str, default=signat_ps.parameters['out_name'].default)
    parser_plot_scatterplot.add_argument('--out_type', type=str, default=signat_ps.parameters['out_type'].default)
    parser_plot_scatterplot.add_argument('--out_directory', type=str, default=signat_ps.parameters['out_directory'].default)
    parser_plot_scatterplot.add_argument('--show', action='store_false', default=signat_ps.parameters['show'].default)
    parser_plot_scatterplot.add_argument('--xlim_kws', type=dict, default=signat_ps.parameters['xlim_kws'].default)
    parser_plot_scatterplot.add_argument('--ylim_kws', type=dict, default=signat_ps.parameters['ylim_kws'].default)
    parser_plot_scatterplot.add_argument('--scatterplot_kws', type=dict, default=signat_ps.parameters['scatterplot_kws'].default)
    parser_plot_scatterplot.add_argument('--subplots_kws', type=dict, default=signat_ps.parameters['subplots_kws'].default)
    parser_plot_scatterplot.add_argument('--axhline_kws', type=dict, default=signat_ps.parameters['axhline_kws'].default)
    parser_plot_scatterplot.set_defaults(func=plot_scatterplot)

    ##################################################

    from be_scan.plot import plot_boxes
    signat_pb = inspect.signature(plot_boxes)
    parser_plot_boxes = subparsers.add_parser('plot_boxes',
                                              help=next(line for line in plot_boxes.__doc__.splitlines() if line),
                                              description=plot_boxes.__doc__,
                                              formatter_class=argparse.RawDescriptionHelpFormatter)
    # required args
    parser_plot_boxes.add_argument('df_filepath', type=str)
    parser_plot_boxes.add_argument('plot_column', type=str)
    parser_plot_boxes.add_argument('-c', '--comparisons', nargs='+', type=str, required=True)
    parser_plot_boxes.add_argument('-pc', '--plot_conditions', nargs='+', type=str, required=True)
    # optional variables
    parser_plot_boxes.add_argument('--filter_val', action='store_true', default=signat_pb.parameters['filter_val'].default)
    parser_plot_boxes.add_argument('--val_cols', nargs='+', type=str, default=signat_pb.parameters['val_cols'].default)
    parser_plot_boxes.add_argument('--val_min', type=float, default=signat_pb.parameters['val_min'].default)
    parser_plot_boxes.add_argument('--filter_params', action='store_true', default=signat_pb.parameters['filter_params'].default)
    parser_plot_boxes.add_argument('--params_cols', nargs='+', type=str, default=signat_pb.parameters['params_cols'].default)
    parser_plot_boxes.add_argument('--params_conditions', nargs='+', type=list, default=signat_pb.parameters['params_conditions'].default)
    parser_plot_boxes.add_argument('--neg_ctrl', action='store_true', default=signat_pb.parameters['neg_ctrl'].default)
    parser_plot_boxes.add_argument('--neg_ctrl_col', type=str, default=signat_pb.parameters['neg_ctrl_col'].default)
    parser_plot_boxes.add_argument('--neg_ctrl_conditions', nargs='+', type=str, default=signat_pb.parameters['neg_ctrl_conditions'].default)
    parser_plot_boxes.add_argument('--xlab', type=str, default=signat_pb.parameters['xlab'].default)
    parser_plot_boxes.add_argument('--ylab', type=str, default=signat_pb.parameters['ylab'].default)
    parser_plot_boxes.add_argument('--col_label', type=str, default=signat_pb.parameters['col_label'].default)
    parser_plot_boxes.add_argument('--savefig', action='store_false', default=signat_pb.parameters['savefig'].default)
    parser_plot_boxes.add_argument('--out_name', type=str, default=signat_pb.parameters['out_name'].default)
    parser_plot_boxes.add_argument('--out_type', type=str, default=signat_pb.parameters['out_type'].default)
    parser_plot_boxes.add_argument('--out_directory', type=str, default=signat_pb.parameters['out_directory'].default)
    parser_plot_boxes.add_argument('--show', action='store_false', default=signat_pb.parameters['show'].default)
    parser_plot_boxes.set_defaults(func=plot_boxes)

    ##################################################

    from be_scan.plot import plot_corr_heatmap
    signat_pch = inspect.signature(plot_corr_heatmap)
    parser_plot_corr_heatmap = subparsers.add_parser('plot_corr_heatmap',
                                                     help=next(line for line in plot_corr_heatmap.__doc__.splitlines() if line),
                                                     description=plot_corr_heatmap.__doc__,
                                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # required args
    parser_plot_corr_heatmap.add_argument('df_filepath', type=str)
    parser_plot_corr_heatmap.add_argument('-c', '--comparisons', nargs='+', type=str, required=True)
    # optional args
    parser_plot_corr_heatmap.add_argument('--corr_type', type=str, default=signat_pch.parameters['corr_type'].default)
    parser_plot_corr_heatmap.add_argument('--filter_val', action='store_true', default=signat_pch.parameters['filter_val'].default)
    parser_plot_corr_heatmap.add_argument('--val_cols', nargs='+', type=str, default=signat_pch.parameters['val_cols'].default)
    parser_plot_corr_heatmap.add_argument('--val_min', type=float, default=signat_pch.parameters['val_min'].default)
    parser_plot_corr_heatmap.add_argument('--filter_params', action='store_true', default=signat_pch.parameters['filter_params'].default)
    parser_plot_corr_heatmap.add_argument('--params_cols', nargs='+', type=str, default=signat_pch.parameters['params_cols'].default)
    parser_plot_corr_heatmap.add_argument('--params_conditions', nargs='+', type=list, default=signat_pch.parameters['params_conditions'].default)
    parser_plot_corr_heatmap.add_argument('--xlab', type=str, default=signat_pch.parameters['xlab'].default)
    parser_plot_corr_heatmap.add_argument('--ylab', type=str, default=signat_pch.parameters['ylab'].default)
    parser_plot_corr_heatmap.add_argument('--title', type=str, default=signat_pch.parameters['title'].default)
    parser_plot_corr_heatmap.add_argument('--out_directory', type=str, default=signat_pch.parameters['out_directory'].default)
    parser_plot_corr_heatmap.add_argument('--out_name', type=str, default=signat_pch.parameters['out_name'].default)
    parser_plot_corr_heatmap.add_argument('--out_type', type=str, default=signat_pch.parameters['out_type'].default)
    parser_plot_corr_heatmap.add_argument('--savefig', action='store_false', default=signat_pch.parameters['savefig'].default)
    parser_plot_corr_heatmap.add_argument('--show', action='store_false', default=signat_pch.parameters['show'].default)
    parser_plot_corr_heatmap.add_argument('--heatmap_kws', type=dict, default=signat_pch.parameters['heatmap_kws'].default)
    parser_plot_corr_heatmap.add_argument('--subplots_kws', type=dict, default=signat_pch.parameters['subplots_kws'].default)
    parser_plot_corr_heatmap.set_defaults(func=plot_corr_heatmap)

    ##################################################

    from be_scan.plot import plot_corr_scatterplot
    signat_pcs = inspect.signature(plot_corr_scatterplot)
    # required args
    parser_plot_corr_scatterplot = subparsers.add_parser('plot_corr_scatterplot',
                                                         help=next(line for line in plot_corr_scatterplot.__doc__.splitlines() if line),
                                                         description=plot_corr_scatterplot.__doc__,
                                                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_plot_corr_scatterplot.add_argument('df_filepath', type=str)
    parser_plot_corr_scatterplot.add_argument('condition1', type=str)
    parser_plot_corr_scatterplot.add_argument('condition2', type=str)
    # optional args
    parser_plot_corr_scatterplot.add_argument('--filter_val', action='store_true', default=signat_pcs.parameters['filter_val'].default)
    parser_plot_corr_scatterplot.add_argument('--val_cols', nargs='+', type=str, default=signat_pcs.parameters['val_cols'].default)
    parser_plot_corr_scatterplot.add_argument('--val_min', type=float, default=signat_pcs.parameters['val_min'].default)
    parser_plot_corr_scatterplot.add_argument('--filter_params', action='store_true', default=signat_pcs.parameters['filter_params'].default)
    parser_plot_corr_scatterplot.add_argument('--params_cols', nargs='+', type=str, default=signat_pcs.parameters['params_cols'].default)
    parser_plot_corr_scatterplot.add_argument('--params_conditions', nargs='+', type=list, default=signat_pcs.parameters['params_conditions'].default)
    parser_plot_corr_scatterplot.add_argument('--include_hue', action='store_true', default=signat_pcs.parameters['include_hue'].default)
    parser_plot_corr_scatterplot.add_argument('--hue_col', type=str, default=signat_pcs.parameters['hue_col'].default)
    parser_plot_corr_scatterplot.add_argument('--hue_order', nargs='+', type=str, default=signat_pcs.parameters['hue_order'].default)
    parser_plot_corr_scatterplot.add_argument('--palette', nargs='+', type=str, default=signat_pcs.parameters['palette'].default)
    parser_plot_corr_scatterplot.add_argument('--savefig', action='store_false', default=signat_pcs.parameters['savefig'].default)
    parser_plot_corr_scatterplot.add_argument('--out_directory', type=str, default=signat_pcs.parameters['out_directory'].default)
    parser_plot_corr_scatterplot.add_argument('--out_name', type=str, default=signat_pcs.parameters['out_name'].default)
    parser_plot_corr_scatterplot.add_argument('--out_type', type=str, default=signat_pcs.parameters['out_type'].default)
    parser_plot_corr_scatterplot.add_argument('--show', action='store_false', default=signat_pcs.parameters['show'].default)
    parser_plot_corr_scatterplot.add_argument('--jointplot_kws', type=dict, default=signat_pcs.parameters['jointplot_kws'].default)
    parser_plot_corr_scatterplot.set_defaults(func=plot_corr_scatterplot)

    ##################################################

### sgrna ###
    
    from be_scan.sgrna import generate_BE_guides
    signat_gBEg = inspect.signature(generate_BE_guides)
    parser_generate_BE_guides = subparsers.add_parser('generate_BE_guides', 
                                                      help=next(line for line in generate_BE_guides.__doc__.splitlines() if line),
                                                      description=generate_BE_guides.__doc__,
                                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_generate_BE_guides.add_argument('gene_filepath', type=str)
    parser_generate_BE_guides.add_argument('cas_type', type=str)
    parser_generate_BE_guides.add_argument('edit_from', type=str)
    parser_generate_BE_guides.add_argument('edit_to', type=str)
    parser_generate_BE_guides.add_argument('--gene_name', type=str, default=signat_gBEg.parameters['gene_name'].default)
    parser_generate_BE_guides.add_argument('--PAM', type=str, default=signat_gBEg.parameters['PAM'].default)
    parser_generate_BE_guides.add_argument('--window', nargs='+', type=int, default=signat_gBEg.parameters['window'].default)
    parser_generate_BE_guides.add_argument('--output_name', type=str, default=signat_gBEg.parameters['output_name'].default)
    parser_generate_BE_guides.add_argument('--output_dir', type=str, default=signat_gBEg.parameters['output_dir'].default)
    parser_generate_BE_guides.add_argument('--return_df', action='store_false', default=signat_gBEg.parameters['return_df'].default)
    parser_generate_BE_guides.add_argument('--save_df', action='store_false', default=signat_gBEg.parameters['save_df'].default)
    parser_generate_BE_guides.add_argument('--exclude_introns', action='store_false', default=signat_gBEg.parameters['exclude_introns'].default)
    parser_generate_BE_guides.add_argument('--exclude_nontargeting', action='store_false', default=signat_gBEg.parameters['exclude_nontargeting'].default)
    parser_generate_BE_guides.set_defaults(func=generate_BE_guides)

    ##################################################
    
    from be_scan.sgrna import check_guides
    signat_cg = inspect.signature(check_guides)
    parser_check_guides = subparsers.add_parser('check_guides', 
                                                help=next(line for line in check_guides.__doc__.splitlines() if line),
                                                description=check_guides.__doc__,
                                                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_check_guides.add_argument('guides_file', type=str)
    parser_check_guides.add_argument('genome_file', type=str)
    parser_check_guides.add_argument('--output_name', type=str, default=signat_cg.parameters['output_name'].default)
    parser_check_guides.add_argument('--output_dir', type=str, default=signat_cg.parameters['output_dir'].default)
    parser_check_guides.add_argument('--delete', action='store_true', default=signat_cg.parameters['delete'].default)
    parser_check_guides.add_argument('--return_df', action='store_false', default=signat_cg.parameters['return_df'].default)
    parser_check_guides.add_argument('--save_df', action='store_false', default=signat_cg.parameters['save_df'].default)
    parser_check_guides.set_defaults(func=check_guides)

    ##################################################
    
    from be_scan.sgrna import annotate_guides
    signat_ag = inspect.signature(annotate_guides)
    parser_annot_guides = subparsers.add_parser('annotate_guides', 
                                                   help=next(line for line in annotate_guides.__doc__.splitlines() if line),
                                                   description=annotate_guides.__doc__,
                                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_annot_guides.add_argument('guides_file', type=str)
    parser_annot_guides.add_argument('edit_from', type=str)
    parser_annot_guides.add_argument('edit_to', type=str)
    parser_annot_guides.add_argument('--gene_filepath', type=str, default=signat_ag.parameters['gene_filepath'].default)
    parser_annot_guides.add_argument('--protein_filepath', type=str, default=signat_ag.parameters['protein_filepath'].default)
    parser_annot_guides.add_argument('--window', nargs='+', type=int, default=signat_ag.parameters['window'].default)
    parser_annot_guides.add_argument('--seq_col', type=str, default=signat_ag.parameters['seq_col'].default)
    parser_annot_guides.add_argument('--gene_pos_col', type=str, default=signat_ag.parameters['gene_pos_col'].default)
    parser_annot_guides.add_argument('--frame_col', type=str, default=signat_ag.parameters['frame_col'].default)
    parser_annot_guides.add_argument('--strand_col', type=str, default=signat_ag.parameters['strand_col'].default)
    parser_annot_guides.add_argument('--output_name', type=str, default=signat_ag.parameters['output_name'].default)
    parser_annot_guides.add_argument('--output_dir', type=str, default=signat_ag.parameters['output_dir'].default)
    parser_annot_guides.add_argument('--return_df', action='store_false', default=signat_ag.parameters['return_df'].default)
    parser_annot_guides.add_argument('--save_df', action='store_false', default=signat_ag.parameters['save_df'].default)
    parser_annot_guides.set_defaults(func=annotate_guides)

    ##################################################
    
    from be_scan.sgrna import guides
    signat_g = inspect.signature(guides)
    parser_guides = subparsers.add_parser('guides', 
                                          help=next(line for line in guides.__doc__.splitlines() if line),
                                          description=guides.__doc__,
                                          formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_guides.add_argument('gene_filepath', type=str)
    parser_guides.add_argument('genome_file', type=str)
    parser_guides.add_argument('cas_type', type=str)
    parser_guides.add_argument('edit_from', type=str)
    parser_guides.add_argument('edit_to', type=str)
    parser_guides.add_argument('--gene_name', type=str, default=signat_g.parameters['gene_name'].default)
    parser_guides.add_argument('--protein_filepath', type=str, default=signat_g.parameters['protein_filepath'].default)
    parser_guides.add_argument('--PAM', type=str, default=signat_g.parameters['PAM'].default)
    parser_guides.add_argument('--window', nargs='+', type=int, default=signat_g.parameters['window'].default)
    parser_guides.add_argument('--output_name', type=str, default=signat_g.parameters['output_name'].default)
    parser_guides.add_argument('--output_dir', type=str, default=signat_g.parameters['output_dir'].default)
    parser_guides.add_argument('--return_df', action='store_false', default=signat_g.parameters['return_df'].default)
    parser_guides.add_argument('--delete', action='store_true', default=signat_g.parameters['delete'].default)
    parser_guides.add_argument('--save_df', action='store_false', default=signat_g.parameters['save_df'].default)
    parser_guides.add_argument('--exclude_introns', action='store_false', default=signat_g.parameters['exclude_introns'].default)
    parser_guides.add_argument('--exclude_nontargeting', action='store_false', default=signat_g.parameters['exclude_nontargeting'].default)
    parser_guides.set_defaults(func=guides)

    args = parser.parse_args()
    print(args)
    function = args.func
    function_args = vars(args)
    del function_args['func']
    print(function_args)
    function(**function_args)
