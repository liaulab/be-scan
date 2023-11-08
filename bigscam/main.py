import argparse
import inspect

def main():
    parser = argparse.ArgumentParser(prog="bigscam")
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

    from .sgrna.findall_be import add_parser_args, main as findall_be_main
    parser_findall_be = subparsers.add_parser("findall_be", description='find all guides accessible for base editing')
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
