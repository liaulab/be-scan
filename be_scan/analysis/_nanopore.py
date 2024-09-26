import re
import itertools
import os
from tqdm import tqdm
import mappy
import os
import multiprocessing as mp
import itertools
from Bio.Seq import reverse_complement
import matplotlib.pyplot as plt
import Bio.SeqIO
import numpy as np
import matplotlib.ticker

def align_helper(args):
    reads_batch, ref_fasta_fname, edit_from_keep, edit_to_keep = args

    aligner = mappy.Aligner(ref_fasta_fname, preset="map-ont")

    out = ""

    for read_id, read_seq, _ in reads_batch:

        # Align the current read
        alignments = aligner.map(read_seq, cs=True)
        alignment = next(alignments, None)

        # Remove multi-mappers and unmapped reads
        if not alignment or next(alignments, None):
            continue

        cur_pos = alignment.r_st
        edit_pos = []
        for cs_op, cs_arg in re.findall(r"([=:*+-~])([\dATCGNatcgn]+)", alignment.cs):
            if cs_op == "*":
                edit_from, edit_to = cs_arg
                edit_from, edit_to = edit_from.upper(), edit_to.upper()
                if edit_from == edit_from_keep and edit_to == edit_to_keep:
                    edit_pos.append(cur_pos)
            cur_pos += ref_bp_consumed(cs_op, cs_arg)

        out += f"{read_id}\t{alignment.ctg}\t{alignment.r_st}\t{alignment.r_en}\t{','.join(map(str, edit_pos))}\n"

    return out

def get_fastq_iter(fastq_file: str) -> tuple[str, str, str]:
    """
    Returns an iterator over the fastq file or directory of fastq files. Each iteration returns read ID, sequence, and quality.
    """
    if os.path.isfile(fastq_file):
        return mappy.fastx_read(fastq_file)
    elif os.path.isdir(fastq_file):
        return itertools.chain.from_iterable(mappy.fastx_read(os.path.join(fastq_file, f)) for f in os.listdir(fastq_file))
    else:
        raise FileNotFoundError(f"fastq file {fastq_file} not found")

def ref_bp_consumed(cs_op: str, cs_arg: str) -> int:
    '''
    Returns the number of bases consumed by the cigar string operation.
    '''
    if cs_op == ":":
        return int(cs_arg)
    elif cs_op == "-":
        return len(cs_arg)
    elif cs_op == "+":
        return 0
    elif cs_op == "*":
        return 1
    raise ValueError("Invalid cigar string")

def find_cutsite(ref_seq: str, guide: str) -> int | None:
    '''
    Zero-indexed cutsite position in ref_seq for guide. Cut is immediately before the base at the indicated position. Assume there's at most one cutsite. If guide is not found, return None.
    '''
    assert len(guide) == 20, "Guide should be 20bp long"
    # + strand
    if (cutsite := ref_seq.find(guide)) != -1:
        return cutsite + 17
    # - strand
    elif (cutsite := ref_seq.find(reverse_complement(guide))) != -1:
        return cutsite + 3
    return None

def align_reads(output_tsv: str, ref_fasta: str, input_fastq: str, edit_from: str, edit_to: str, plot_fname: str | None = None, guides: str | None = None):
    """
    Uses minimap to align fastq reads to reference fasta. Output TSV file contains read ID, name of reference sequence, start and end positions of alignment, and a comma-separated list of positions of edits. Positions are zero-indexed. Start position is inclusive, end position is exclusive.
    :param output_tsv: Path to output TSV file
    :param ref_fasta: Path to reference fasta file
    :param input_fastq: Path to input fastq file. Or path to directory containing fastq files.
    :param edit_from: Check for edits from this base
    :param edit_to: Check for edits to this base
    :param plot_fname: Path to save plot. If None, no plot is saved.
    :param guides: Comma-separated list of guides to plot. If None, no guides are plotted.
    """

    edit_from, edit_to = edit_from.upper(), edit_to.upper()

    with open(output_tsv, "wt") as aln_file, mp.Pool(mp.cpu_count()) as pool:
        aln_file.write("read_id\tref_name\tstart\tend\tedit_pos\n")
        for out in pool.imap(
                align_helper,
                ((reads_batch, ref_fasta, edit_from, edit_to) for reads_batch in itertools.batched(tqdm(get_fastq_iter(input_fastq)), 10000)),
                ):
            if out is not None:
                aln_file.write(out)

    if plot_fname is not None:
        # Read in reference sequences
        ref_seqs = {record.id: str(record.seq) for record in Bio.SeqIO.parse(ref_fasta, "fasta")}
        # Count edits per position
        edit_pos_hist = {name: np.zeros(len(seq), dtype=int) for name, seq in ref_seqs.items()}
        with open(output_tsv, "rt") as aln_file:
            next(aln_file, None) # skip header
            for line in aln_file:
                _, ref_name, _, _, edit_pos = line.rstrip("\n").split("\t")
                if edit_pos == "":
                    continue
                edits = map(int, edit_pos.split(","))
                for edit in edits:
                    edit_pos_hist[ref_name][edit] += 1
        cutsites = {name: [] for name in ref_seqs}
        if guides is not None:
            # Read guides
            for guide in guides.split(","):
                for name, seq in ref_seqs.items():
                    if (cutsite := find_cutsite(seq, guide)) is not None:
                        cutsites[name].append(cutsite)

        fig, axes = plt.subplots(figsize=(12, 4 * len(edit_pos_hist)))
        if len(edit_pos_hist) == 1:
            axes = [axes]
        for ax, (name, _edit_pos_hist) in zip(axes, edit_pos_hist.items()):
            ax.set_title(name)
            ax.plot(_edit_pos_hist)
            for cut_site in cutsites[name]:
                ax.axvline(x=cut_site, color="red", linestyle="dashed")
            ax.set_ylabel("Edit coverage")
            ax.set_xlabel("Position in the locus")
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=max(_edit_pos_hist.sum(), 1)))
        fig.savefig(plot_fname, bbox_inches="tight", dpi=300)
        plt.close(fig)
