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
    parser_count_reads.add_argument('--DIR', choices=["FWD", "REV"], default=signature_count_reads.parameters['DIR'].default)
    parser_count_reads.add_argument('--KEY', type=str, default=signature_count_reads.parameters['KEY'].default)
    parser_count_reads.add_argument('--KEY_REV', type=str, default=signature_count_reads.parameters['KEY_REV'].default)
    parser_count_reads.add_argument('--out_counts', type=str, default=signature_count_reads.parameters['out_counts'].default)
    parser_count_reads.add_argument('--out_np', type=str, default=signature_count_reads.parameters['out_np'].default)
    parser_count_reads.add_argument('--out_stats', type=str, default=signature_count_reads.parameters['out_stats'].default)
    parser_count_reads.set_defaults(func=count_reads)

    args = parser.parse_args()
    function = args.func
    function_args = vars(args)
    del function_args['func']
    function(**function_args)
