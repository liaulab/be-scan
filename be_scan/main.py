import argparse
import inspect

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def main():

### analysis ###

    parser = argparse.ArgumentParser(prog="be_scan")
    subparsers = parser.add_subparsers(required=True)

    from .analysis import count_reads
    signature_count_reads = inspect.signature(count_reads)
    parser_count_reads = subparsers.add_parser('count_reads',
                                               help=next(line for line in count_reads.__doc__.splitlines() if line),
                                               description=count_reads.__doc__,
                                               formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_count_reads.add_argument('in_fastq', type=str)
    parser_count_reads.add_argument('in_ref', type=str)
    parser_count_reads.add_argument('--KEY_INTERVAL', type=int, nargs=2, default=signature_count_reads.parameters['KEY_INTERVAL'].default)
    parser_count_reads.add_argument('--KEY', type=str, default=signature_count_reads.parameters['KEY'].default)
    parser_count_reads.add_argument('--KEY_REV', type=str, default=signature_count_reads.parameters['KEY_REV'].default)
    parser_count_reads.add_argument('--out_counts', type=str, default=signature_count_reads.parameters['out_counts'].default)
    parser_count_reads.add_argument('--out_np', type=str, default=signature_count_reads.parameters['out_np'].default)
    parser_count_reads.add_argument('--out_stats', type=str, default=signature_count_reads.parameters['out_stats'].default)
    parser_count_reads.add_argument('--dont_trim_G', action='store_true')
    parser_count_reads.set_defaults(func=count_reads)

    ##################################################

    from .analysis import validate_cloning
    signature_validate_cloning = inspect.signature(validate_cloning)
    parser_validate_cloning = subparsers.add_parser('validate_cloning',
                                                    help=next(line for line in validate_cloning.__doc__.splitlines() if line),
                                                    description=validate_cloning.__doc__,
                                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser_validate_cloning.add_argument('query_dir', type=str)
    parser_validate_cloning.add_argument('spacers', type=str)
    parser_validate_cloning.add_argument('vector', type=str)
    parser_validate_cloning.add_argument('out_csv', type=str, help="Output CSV file name")
    parser_validate_cloning.add_argument('--enzyme', type=str, default="Esp3I", help="Default: %(default)s")
    parser_validate_cloning.add_argument('--overhangs', type=str, default=('CACC', 'GTTT'), nargs=2, help="Default: %(default)s")
    parser_validate_cloning.add_argument('--ref_name_pattern', type=str, default=signature_validate_cloning.parameters['ref_name_pattern'].default, help="Default: %(default)s")
    parser_validate_cloning.add_argument('--flank_width', type=int, default=signature_validate_cloning.parameters['flank_width'].default, help="Default: %(default)s")
    parser_validate_cloning.set_defaults(func=validate_cloning)

    ##################################################

### plot ###

    from .plot.scatterplot import plot_scatterplot
    # signature_plot_scatterplot = inspect.signature(plot_scatterplot)
    parser_plot_scatterplot = subparsers.add_parser('plot_scatterplot', 
                                                    description='This function takes in a dataframe from count_reads, performs normalization, and then plots the data for each condition to reveal which guides are enriched')
    parser_plot_scatterplot.add_argument('-df', '--df_filepath', type=str, help='filepath to .csv data generated from count_reads')
    parser_plot_scatterplot.add_argument('-x', '--x_column', type=str, help='column of .csv for x axis, typically amino acid position')
    parser_plot_scatterplot.add_argument('-y', '--y_column', type=str, help='column of .csv for y axis, typically the normalized log_fc change score')
    parser_plot_scatterplot.add_argument('-hue', '--hue_column', type=str, help='column of .csv which correspond to coloring of scatterplot points')
    parser_plot_scatterplot.add_argument('-ncat', '--neg_ctrl_col', type=str, help='column of .csv which correspond to normalization control')
    parser_plot_scatterplot.add_argument('-ncol', '--neg_ctrl_category', type=str, help='categorical variable of neg_ctrl_col .csv which correspond to normalization control')
    parser_plot_scatterplot.add_argument('-c','--comparisons', nargs='+', type=str, help='list of comparisons that correspond to columns of .csv data', required=True)
    
    # optional variables
    parser_plot_scatterplot.add_argument('-win','--window', nargs='+', type=int, help='inclusive window of which amino acid positions are shown in the plot', default=None)
    
    # output variables (all optional: default is './scatterplot.png')
    parser_plot_scatterplot.add_argument('--x_label', type=str, default='Amino Acid Position', help="name of the x-axis label")
    parser_plot_scatterplot.add_argument('--y_label', type=str, default='sgRNA Score', help="name of y-axis label")
    parser_plot_scatterplot.add_argument('--out_name', type=str, default='scatterplot', help="name of figure output")
    parser_plot_scatterplot.add_argument('--out_type', type=str, default='pdf', help='file type of figure output')
    parser_plot_scatterplot.add_argument('--out_directory', type=str, default='', help="path to output directory")
    parser_plot_scatterplot.add_argument('--alpha', type=float, default=0.8, help="transparency of scatterplot points")
    parser_plot_scatterplot.add_argument('--linewidth', type=float, default=1.0, help="linewidth of plot")
    parser_plot_scatterplot.add_argument('--edgecolor', type=str, default='black', help="color of scatterplot edge lines")
    parser_plot_scatterplot.add_argument('--s', type=float, default=25, help="size of scatterplot points")
    parser_plot_scatterplot.add_argument('--figsize', nargs='+', type=int, help='the figsize (length, width)', default=[8, 4])
    parser_plot_scatterplot.set_defaults(func=plot_scatterplot)

    ##################################################

    from .plot.correlation_heatmap import plot_corr_heatmap
    # signature_plot_corr_heatmap = inspect.signature(plot_corr_heatmap)
    parser_plot_corr_heatmap = subparsers.add_parser('plot_corr_heatmap', 
                                                     description='his function takes in a dataframe from count_reads, and plots a scatterplot showing correlation between two given conditions')
    parser_plot_corr_heatmap.add_argument('-df', '--df_filepath', type=str, help='filepath to .csv data generated from count_reads')
    parser_plot_corr_heatmap.add_argument('-c','--comparisons', nargs='+', type=str, help='list of comparisons that correspond to columns of .csv data', required=True)
    
    # output variables (all optional: default is './correlation_heatmap.png')
    parser_plot_corr_heatmap.add_argument('-ct', '--corr_type', type=str, default='spearman', help='type of correlation calculation')
    parser_plot_corr_heatmap.add_argument('--xlab', type=str, default='', help="name of the x-axis label")
    parser_plot_corr_heatmap.add_argument('--ylab', type=str, default='', help="name of y-axis label")
    parser_plot_corr_heatmap.add_argument('--title', type=str, default='', help="name of title label")
    parser_plot_corr_heatmap.add_argument('--out_name', type=str, default='correlation_heatmap', help="name of figure output")
    parser_plot_corr_heatmap.add_argument('--out_type', type=str, default='pdf', help='file type of figure output')
    parser_plot_corr_heatmap.add_argument('--out_directory', type=str, default='', help="path to output directory")
    parser_plot_corr_heatmap.add_argument('--center', type=float, default=0, help="value for centering the colormap")
    parser_plot_corr_heatmap.add_argument('--linewidth', type=float, default=0.5, help="linewidth of lines dividing cells")
    parser_plot_corr_heatmap.add_argument('--cmap', type=str, default='coolwarm', help="positions to put larger lines for dividing the cells")
    parser_plot_corr_heatmap.add_argument('--line_pos', nargs='+', type=int, help='positions to put larger lines for dividing the cells')
    parser_plot_corr_heatmap.set_defaults(func=plot_corr_heatmap)

    ##################################################

    from .plot.correlation_scatter import plot_corr_scatterplot
    # signature_plot_corr_scatterplot = inspect.signature(plot_corr_scatterplot)
    parser_plot_corr_scatterplot = subparsers.add_parser('plot_corr_scatterplot', 
                                                    description='This function takes in a dataframe from count_reads, and plots a heatmap showing correlation between all given comparison conditions')
    parser_plot_corr_scatterplot.add_argument('-df', '--df_filepath', type=str, help='filepath to .csv data generated from count_reads')
    parser_plot_corr_scatterplot.add_argument('-c1', '--condition1', type=str, help='comparison condition 1, name of a column in .csv data')
    parser_plot_corr_scatterplot.add_argument('-c2', '--condition2', type=str, help='comparison condition 2, name of a column in .csv data')
    parser_plot_corr_scatterplot.add_argument('-hue', '--hue_column', type=str, help='the categorial data for each of the points, name of a column in .csv data')
    parser_plot_corr_scatterplot.add_argument('--hue_order', nargs='+', type=str, help='a list of categorial variables in hue_column')
    parser_plot_corr_scatterplot.add_argument('--palette', nargs='+', type=str, help='a list of colors which correspond to hue_order')
    
    parser_plot_corr_scatterplot.add_argument('--xmin', type=float, default=None, help="x-axis left bound")
    parser_plot_corr_scatterplot.add_argument('--xmax', type=float, default=None, help="x-axis right bound")
    parser_plot_corr_scatterplot.add_argument('--ymin', type=float, default=None, help="y-axis lower bound")
    parser_plot_corr_scatterplot.add_argument('--ymax', type=float, default=None, help="y-axis upper bound")
    parser_plot_corr_scatterplot.add_argument('--xlab', type=str, default='cond1 score', help="name of the x-axis label")
    parser_plot_corr_scatterplot.add_argument('--ylab', type=str, default='cond2 score', help="name of y-axis label")
    # output variables (all optional: default is './correlation_scatterplot.png')
    parser_plot_corr_scatterplot.add_argument('--out_name', type=str, default='correlation_scatterplot', help="name of figure output")
    parser_plot_corr_scatterplot.add_argument('--out_type', type=str, default='pdf', help='file type of figure output')
    parser_plot_corr_scatterplot.add_argument('--out_directory', type=str, default='', help="path to output directory")
    parser_plot_corr_scatterplot.add_argument('--alpha', type=float, default=0.8, help="transparency of scatterplot points")
    parser_plot_corr_scatterplot.add_argument('--linewidth', type=float, default=1.0, help="linewidth of plot")
    parser_plot_corr_scatterplot.add_argument('--edgecolor', type=str, default='black', help="color of scatterplot edge lines")
    parser_plot_corr_scatterplot.add_argument('--s', type=float, default=25, help="size of scatterplot points")
    parser_plot_corr_scatterplot.add_argument('--figsize', nargs='+', type=int, help='the figsize (length, width)', default=[4.5, 4])
    parser_plot_corr_scatterplot.set_defaults(func=plot_corr_scatterplot)

    ##################################################



    ##################################################

### sgrna ###

    from .sgrna.findall_be import add_parser_args, main as findall_be_main
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
