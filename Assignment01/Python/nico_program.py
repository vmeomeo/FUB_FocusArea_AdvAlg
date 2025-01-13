import gzip
import iv2py as iv
directory='/home/nico/Desktop/University/1AdvancedAlgo/assignment1/FUB_FocusArea_AdvAlg/'
hg=gzip.open(directory+'/Assignment01/data/hg38_partial.fasta.gz','r')
read_40=gzip.open(directory+'/Assignment01/data/illumina_reads_40.fasta.gz','r')


for record in iv.fasta.reader(file=directory+'/Assignment01/data/hg38_partial.fasta.gz'):
    sequence= record.seq

reads=[]
for record in iv.fasta.reader(file=directory+'/Assignment01/data/illumina_reads_40.fasta.gz'):
    reads.append(record.seq)


import time

# Start the timer
start_time = time.time()
for i in range(10):
    pos = []
    start = 0
    seq = reads[i]
    while True:
        start = sequence.find(seq, start)
        if start == -1:  # No more matches
            break
        pos.append(start)
        start += 1

end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)


