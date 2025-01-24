import iv2py as iv

def read_fasta(file_path):
    sequence = []
    for record in iv.fasta.reader(file_path):
        sequence.append(record.seq)  # Extract only the sequence part
    return sequence
reference= read_fasta("../data/hg38_partial.fasta.gz")
#print(len(reference[0])) 100M


def create_suffix_array(reference_file):
    sa= iv.create_suffixarray(reference[0])
    return sa
suf_arr= create_suffix_array(reference)
