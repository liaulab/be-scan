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
    parser_count_reads.add_argument('in_ref', type=str)
    parser_count_reads.add_argument('--file_dir', type=str, default=signat_cr.parameters['file_dir'].default)
    parser_count_reads.add_argument('--KEY_INTERVAL', type=int, nargs=2, default=signat_cr.parameters['KEY_INTERVAL'].default)
    parser_count_reads.add_argument('--KEY', type=str, default=signat_cr.parameters['KEY'].default)
    parser_count_reads.add_argument('--KEY_REV', type=str, default=signat_cr.parameters['KEY_REV'].default)
    parser_count_reads.add_argument('--dont_trim_G', action='store_true')
    parser_count_reads.set_defaults(func=count_reads)

    ##################################################
    
    from be_scan.analysis.merge_and_norm import merge_and_norm
    signat_man = inspect.signature(merge_and_norm)
    parser_merge_and_norm = subparsers.add_parser('merge_and_norm', 
                                                  help=next(line for line in merge_and_norm.__doc__.splitlines() if line),
                                                  description=merge_and_norm.__doc__,
                                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_merge_and_norm.add_argument('sample_sheet', type=str)
    parser_merge_and_norm.add_argument('in_ref', type=str)
    parser_merge_and_norm.add_argument('--file_dir', type=str, default=signat_man.parameters['file_dir'].default)
    parser_merge_and_norm.add_argument('--t0', type=str, default=signat_man.parameters['t0'].default)
    parser_merge_and_norm.add_argument('--dir_counts', type=str, default=signat_man.parameters['dir_counts'].default)
    parser_merge_and_norm.add_argument('--save', type=str, default=signat_man.parameters['save'].default)
    parser_merge_and_norm.add_argument('--out_reads', type=str, default=signat_man.parameters['out_reads'].default)
    parser_merge_and_norm.add_argument('--out_log2', type=str, default=signat_man.parameters['out_log2'].default)
    parser_merge_and_norm.add_argument('--out_t0', type=str, default=signat_man.parameters['out_t0'].default)
    parser_merge_and_norm.add_argument('--return_df', type=bool, default=signat_man.parameters['return_df'].default)
    parser_merge_and_norm.set_defaults(func=merge_and_norm)
    
    ##################################################
    
    from be_scan.analysis.average_reps import average_reps
    signat_ar = inspect.signature(average_reps)
    parser_average_reps = subparsers.add_parser('average_reps', 
                                                help=next(line for line in average_reps.__doc__.splitlines() if line),
                                                description=average_reps.__doc__,
                                                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_average_reps.add_argument('sample_sheet', type=str)
    parser_average_reps.add_argument('in_lfc', type=str)
    parser_average_reps.add_argument('--file_dir', type=str, default=signat_ar.parameters['file_dir'].default)
    parser_average_reps.add_argument('--save', type=bool, default=signat_ar.parameters['save'].default)
    parser_average_reps.add_argument('--out_conds', type=str, default=signat_ar.parameters['out_conds'].default)
    parser_average_reps.add_argument('--return_df', type=bool, default=signat_ar.parameters['return_df'].default)
    parser_average_reps.set_defaults(func=average_reps)
    
    ##################################################

### plot ###

    from be_scan.plot.scatterplot import plot_scatterplot
    signat_ps = inspect.signature(plot_scatterplot)
    parser_plot_scatterplot = subparsers.add_parser('plot_scatterplot', 
                                                    help=next(line for line in plot_scatterplot.__doc__.splitlines() if line),
                                                    description=plot_scatterplot.__doc__,
                                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_plot_scatterplot.add_argument('-df', '--df_filepath', type=str, required=True)
    parser_plot_scatterplot.add_argument('-x', '--x_column', type=str, required=True)
    parser_plot_scatterplot.add_argument('-y', '--y_column', type=str, required=True)
    parser_plot_scatterplot.add_argument('-hue', '--hue_column', type=str, required=True)
    parser_plot_scatterplot.add_argument('-c','--comparisons', nargs='+', type=str, required=True)
    parser_plot_scatterplot.add_argument('-ncol', '--neg_ctrl_col', type=str, required=True)
    parser_plot_scatterplot.add_argument('-ncat', '--neg_ctrl_category', type=str, required=True)
    # optional variables
    parser_plot_scatterplot.add_argument('--xmin', type=float, default=signat_ps.parameters['xmin'].default)
    parser_plot_scatterplot.add_argument('--xmax', type=float, default=signat_ps.parameters['xmax'].default)    
    parser_plot_scatterplot.add_argument('--xlab', type=str, default=signat_ps.parameters['xlab'].default)
    parser_plot_scatterplot.add_argument('--ylab', type=str, default=signat_ps.parameters['ylab'].default)
    parser_plot_scatterplot.add_argument('--out_name', type=str, default=signat_ps.parameters['out_name'].default)
    parser_plot_scatterplot.add_argument('--out_type', type=str, default=signat_ps.parameters['out_type'].default)
    parser_plot_scatterplot.add_argument('--out_directory', type=str, default=signat_ps.parameters['out_directory'].default)
    parser_plot_scatterplot.add_argument('--savefig', type=bool, default=signat_ps.parameters['savefig'].default)
    parser_plot_scatterplot.set_defaults(func=plot_scatterplot)

    ##################################################

    from be_scan.plot.boxes import plot_boxes
    signat_pb = inspect.signature(plot_boxes)
    parser_plot_boxes = subparsers.add_parser('plot_boxes',
                                              help=next(line for line in plot_boxes.__doc__.splitlines() if line),
                                              description=plot_boxes.__doc__,
                                              formatter_class=argparse.RawDescriptionHelpFormatter)
    # required args
    parser_plot_boxes.add_argument('-df', '--df_filepath', type=str, required=True)
    parser_plot_boxes.add_argument('-p', '--plot_column', type=str, required=True)
    parser_plot_boxes.add_argument('-pc', '--plot_conditions', nargs='+', type=str, required=True)
    parser_plot_boxes.add_argument('-y', '--y_column', type=str, required=True)
    parser_plot_boxes.add_argument('-c', '--comparisons', nargs='+', type=str, required=True)
    parser_plot_boxes.add_argument('-ncol', '--neg_ctrl_col', type=str, required=True)
    parser_plot_boxes.add_argument('-ncat', '--neg_ctrl_category', type=str, required=True)
    # optional variables
    parser_plot_boxes.add_argument('--filter_column', type=str, default=signat_pb.parameters['filter_column'].default)
    parser_plot_boxes.add_argument('--filter_category', type=str, default=signat_pb.parameters['filter_category'].default)
    parser_plot_boxes.add_argument('--xlab', type=str, default=signat_pb.parameters['xlab'].default)
    parser_plot_boxes.add_argument('--ylab', type=str, default=signat_pb.parameters['ylab'].default)
    parser_plot_boxes.add_argument('--out_name', type=str, default=signat_pb.parameters['out_name'].default)
    parser_plot_boxes.add_argument('--out_type', type=str, default=signat_pb.parameters['out_type'].default)
    parser_plot_boxes.add_argument('--out_directory', type=str, default=signat_pb.parameters['out_directory'].default)
    parser_plot_boxes.add_argument('--savefig', type=bool, default=signat_pb.parameters['savefig'].default)
    ### extra args
    parser_plot_boxes.set_defaults(func=plot_boxes)

    ##################################################

    from be_scan.plot.correlation_heatmap import plot_corr_heatmap
    signat_pch = inspect.signature(plot_corr_heatmap)
    parser_plot_corr_heatmap = subparsers.add_parser('plot_corr_heatmap',
                                                     help=next(line for line in plot_corr_heatmap.__doc__.splitlines() if line),
                                                     description=plot_corr_heatmap.__doc__,
                                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # required args
    parser_plot_corr_heatmap.add_argument('-df', '--df_filepath', type=str, required=True)
    parser_plot_corr_heatmap.add_argument('-c','--comparisons', nargs='+', type=str, required=True)
    # optional args
    parser_plot_corr_heatmap.add_argument('--corr_type', type=str, default=signat_pch.parameters['corr_type'].default)
    parser_plot_corr_heatmap.add_argument('--xlab', type=str, default=signat_pch.parameters['xlab'].default)
    parser_plot_corr_heatmap.add_argument('--ylab', type=str, default=signat_pch.parameters['ylab'].default)
    parser_plot_corr_heatmap.add_argument('--title', type=str, default=signat_pch.parameters['title'].default)
    parser_plot_corr_heatmap.add_argument('--out_directory', type=str, default=signat_pch.parameters['out_directory'].default)
    parser_plot_corr_heatmap.add_argument('--out_name', type=str, default=signat_pch.parameters['out_name'].default)
    parser_plot_corr_heatmap.add_argument('--out_type', type=str, default=signat_pch.parameters['out_type'].default)
    parser_plot_corr_heatmap.add_argument('--savefig', type=bool, default=signat_pch.parameters['savefig'].default)
    parser_plot_corr_heatmap.set_defaults(func=plot_corr_heatmap)

    ##################################################

    from be_scan.plot.correlation_scatter import plot_corr_scatterplot
    from be_scan.plot._annotating_ import color_list, list_muttypes
    signat_pcs = inspect.signature(plot_corr_scatterplot)
    # required args
    parser_plot_corr_scatterplot = subparsers.add_parser('plot_corr_scatterplot',
                                                         help=next(line for line in plot_corr_scatterplot.__doc__.splitlines() if line),
                                                         description=plot_corr_scatterplot.__doc__,
                                                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_plot_corr_scatterplot.add_argument('-df', '--df_filepath', type=str, required=True)
    parser_plot_corr_scatterplot.add_argument('-c1', '--condition1', type=str, required=True)
    parser_plot_corr_scatterplot.add_argument('-c2', '--condition2', type=str, required=True)
    parser_plot_corr_scatterplot.add_argument('-hue', '--hue_column', type=str, required=True)
    # optional args
    parser_plot_corr_scatterplot.add_argument('--hue_order', nargs='+', type=str, default=signat_pcs.parameters['hue_order'].default)
    parser_plot_corr_scatterplot.add_argument('--palette', nargs='+', type=str, default=signat_pcs.parameters['palette'].default)
    parser_plot_corr_scatterplot.add_argument('--xmin', type=float, default=signat_pcs.parameters['xmin'].default)
    parser_plot_corr_scatterplot.add_argument('--xmax', type=float, default=signat_pcs.parameters['xmax'].default)
    parser_plot_corr_scatterplot.add_argument('--ymin', type=float, default=signat_pcs.parameters['ymin'].default)
    parser_plot_corr_scatterplot.add_argument('--ymax', type=float, default=signat_pcs.parameters['ymax'].default)
    parser_plot_corr_scatterplot.add_argument('--xlab', type=str, default=signat_pcs.parameters['xlab'].default)
    parser_plot_corr_scatterplot.add_argument('--ylab', type=str, default=signat_pcs.parameters['ylab'].default)
    parser_plot_corr_scatterplot.add_argument('--out_name', type=str, default=signat_pcs.parameters['out_name'].default)
    parser_plot_corr_scatterplot.add_argument('--out_type', type=str, default=signat_pcs.parameters['out_type'].default)
    parser_plot_corr_scatterplot.add_argument('--out_directory', type=str, default=signat_pcs.parameters['out_directory'].default)
    parser_plot_corr_scatterplot.add_argument('--savefig', type=bool, default=signat_pcs.parameters['savefig'].default)
    parser_plot_corr_scatterplot.set_defaults(func=plot_corr_scatterplot)

    ##################################################

### sgrna ###

    from be_scan.sgrna.findall_be import add_parser_args, main as findall_be_main
    parser_findall_be = subparsers.add_parser("findall_be", 
                                              description='find all guides accessible for base editing')
    parser_findall_be = add_parser_args(parser_findall_be)
    parser_findall_be.set_defaults(func=findall_be_main)

    args = parser.parse_args()
    if args.func == findall_be_main:
        findall_be_main(args)
    else:
        function = args.func
        function_args = vars(args)
        del function_args['func']
        function(**function_args)
