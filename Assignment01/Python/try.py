import iv2py as iv

reference = "../data/illumina_reads_80.fasta.gz"

def read_fasta(file_path):
    sequence=[]
    for record in iv.fasta.reader(file_path):
        sequence.append(record)
    return sequence
fasta=read_fasta(reference)
sa = iv.create_suffixarray(fasta)

