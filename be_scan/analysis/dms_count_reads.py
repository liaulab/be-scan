from Bio import SeqIO

def count_spacers_fwd(
    in_lib, in_fastq, column_name,
    nucleotide_column = 'NucleotideSeq',
    START_KEY = "GGTGGAGACTTGTATGTGGTG", # identifies constant seq before mutated region to determine position
    END_KEY = "ATGTGGAAGTGCAACAATGCCACC",  # identifies constant seq after mutated region
    ):

    dict_perfects = {x:0 for x in in_lib['NucleotideSeq']}
    readiter = SeqIO.parse(in_fastq, "fastq")

    for record in readiter: # contains the seq and Qscore etc.
        read_sequence = str.upper(str(record.seq))
        KEY_REGION_START = read_sequence.find(START_KEY)
        KEY_REGION_END = read_sequence.find(END_KEY)
        if (KEY_REGION_START >=0) and (KEY_REGION_END >=0):
            guide = read_sequence[(KEY_REGION_START+len(START_KEY)):(KEY_REGION_END)]
            if guide in dict_perfects.keys():
                dict_perfects[guide] += 1

    in_lib[column_name] = in_lib[nucleotide_column].apply(lambda x: dict_perfects[x])
    return in_lib
