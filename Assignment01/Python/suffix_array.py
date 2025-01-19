# Main logic
import iv

# Read the genome sequence
for record in iv.fasta.reader(file=directory + '/Assignment01/data/hg38_partial.fasta.gz'):
    sequence = record.seq

# Build the suffix array for the genome sequence
suffix_array = build_suffix_array(sequence)

# Read the reads
reads = []
for record in iv.fasta.reader(file=directory + '/Assignment01/data/illumina_reads_40.fasta.gz'):
    reads.append(record.seq)

# Find positions for each read
for i in range(10):  # Limit to the first 10 reads for demonstration
    read = reads[i]
    positions = binary_search_suffix_array(sequence, suffix_array, read)
    print(f"Read {i}: Found at positions {positions}")

def build_suffix_array(sequence):
    """Builds the suffix array for the given sequence."""
    return sorted(range(len(sequence)), key=lambda i: sequence[i:])

def binary_search_suffix_array(sequence, suffix_array, read):
    """Performs binary search to find the range of matches for the read in the suffix array."""
    lo, hi = 0, len(suffix_array)
    # Find the lower bound
    while lo < hi:
        mid = (lo + hi) // 2
        if sequence[suffix_array[mid]:].startswith(read): #If it found already the right sequence
            hi = mid
        else:
            if sequence[suffix_array[mid]:] < read:
                lo = mid + 1
            else:
                hi = mid
    start = lo

    # Find the upper bound
    lo, hi = start, len(suffix_array)
    while lo < hi:
        mid = (lo + hi) // 2
        if sequence[suffix_array[mid]:].startswith(read):
            lo = mid + 1
        else:
            hi = mid
    end = lo
    return suffix_array[start:end]
