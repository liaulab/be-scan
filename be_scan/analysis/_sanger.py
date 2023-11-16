import Bio.SeqIO
import Bio.Seq
import os
import numpy as np
import re
import pandas as pd
import warnings

def validate_cloning(query_dir, spacers, vector, enzyme, overhangs, ref_name_pattern=r"^(.*)[-_]", flank_width=50, out_csv=None):
    """
    Parameters:
        query_dir: directory containing .ab1 files
        spacers: pd.Series with index of plasmid names and values of spacer sequences. Or if string, path to .csv 
                 file with columns "plasmid" and "spacer"
        vector: filename of vector backbone in GenBank format
        enzyme: name of enzyme used to cut the vector. Must be one of "Esp3I", "BpiI", "BbsI", or "BsmBI"
        overhangs: tuple of overhangs used to ligate the spacer into the vector. e.g. ("CACC", "GTTT") for Esp3I
        ref_name_pattern: regular expression to extract plasmid name from filename
        flank_width: number of bases from the backbone to include on either side of the spacer in the reference sequence
    Returns:
        pd.DataFrame with columns "spacer", "reference plasmid", and "success" and index of sanger trace filenames
    """
    flank_left, flank_right = _golden_gate(vector, enzyme, overhangs)
    if isinstance(spacers, str):
        spacers = pd.read_csv(spacers, index_col="plasmid")["spacer"]
    reference_seqs = str(flank_left[-flank_width:]) + spacers + str(flank_right[:flank_width])
    results = []
    for fname in os.listdir(query_dir):
        try:
            sanger_trace = next(Bio.SeqIO.parse(os.path.join(query_dir, fname), "abi"))
            plasmid_name = re.match(ref_name_pattern, fname).group(1)
        except:
            warnings.warn("Could not parse %s" % fname)
            continue
        if plasmid_name not in reference_seqs:
            warnings.warn("No reference sequence for %s" % fname)
            continue
        moving_average_quality = np.convolve(sanger_trace.letter_annotations["phred_quality"], np.ones(10), mode="same") / 10
        trim = np.nonzero(moving_average_quality > 30)[0]
        sanger_trace = sanger_trace[trim.min():trim.max()]
        results.append((fname, plasmid_name, reference_seqs.loc[plasmid_name] in sanger_trace.seq))
    out = pd.merge(spacers.rename("spacer"), pd.DataFrame.from_records(results, columns=["sanger trace fname", "reference plasmid", "success"]), left_index=True, right_on="reference plasmid").set_index("sanger trace fname")
    if out_csv is not None:
        out.to_csv(out_csv)
    return out

def _golden_gate(vector, enzyme, overhangs):
    if enzyme in ("BpiI", "BbsI"):
        cut_seq = Bio.Seq.Seq("GAAGAC")
        cut_site = 2
    elif enzyme in ("Esp3I", "BsmBI"):
        cut_seq = Bio.Seq.Seq("CGTCTC")
        cut_site = 1
    else:
        raise ValueError("Unknown enzyme: %s" % enzyme)

    vector_seq = Bio.SeqIO.read(vector, "genbank").seq

    # plasmids are circular, so put the cut site in the middle of the vector
    while True:
        cutsite_5p = vector_seq.find(cut_seq.reverse_complement()) - cut_site
        cutsite_3p = vector_seq.find(cut_seq) + cut_site + len(cut_seq)
        if cutsite_5p < cutsite_3p:
            break
        vector_seq = vector_seq[1:] + vector_seq[0]
    bp_to_roll = (cutsite_5p + cutsite_3p) // 2 - len(vector_seq) // 2
    vector_seq = vector_seq[bp_to_roll:] + vector_seq[:bp_to_roll]
    cutsite_5p -= bp_to_roll
    cutsite_3p -= bp_to_roll

    flank_left = vector_seq[:cutsite_5p]
    flank_right = vector_seq[cutsite_3p:]
    if flank_left.endswith(overhangs[0]) and flank_right.startswith(overhangs[1]):
        return flank_left, flank_right
    else:
        return flank_right.reverse_complement(), flank_left.reverse_complement()
