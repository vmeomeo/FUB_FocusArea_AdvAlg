import argparse
import gzip
import iv2py as iv

def read_fasta(file_path):
    sequence=[]
    for record in iv.fasta.reader(file_path):
        sequence.append(record)
    return sequence




def binary_search_suffix_array(sequence, suffix_array, read):
    """Performs binary search to find the range of matches for the read in the suffix array."""
    lo, hi = 0, len(suffix_array)
    # Find the lower bound
    while lo < hi:
        mid = (lo + hi) // 2
        if sequence[suffix_array[mid]:].startswith(read):
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

def main():
    parser = argparse.ArgumentParser(description="Search algorithm for DNA sequences using a portion of the reference.")
    parser.add_argument("-reference", type=str, required=True,
                        help="Path to the file containing the reference sequence (FASTA format).")
    parser.add_argument("-query", type=str, required=True,
                        help="Path to the file containing the query sequences (FASTA format).")
    parser.add_argument("-query_ct", type=int, required=True,
                        help="Number of query sequences to search for.")

    args = parser.parse_args()

    try:
        # Read reference and query sequences
        reference_sequence = read_fasta(args.reference)
        query_sequences = read_fasta(args.query)

        # Combine all reference sequences into a single string

        sa = iv.create_suffixarray(reference_sequence)
        print('Suffix array done')


        # Ensure query count does not exceed available queries
        if args.query_ct > len(query_sequences):
            print(f"Error: The specified query count ({args.query_ct}) exceeds the number of queries available ({len(query_sequences)}).")
            return

        # Perform search for each query sequence
        for i in range(args.query_ct):
            query = query_sequences[i]
            positions = binary_search_suffix_array(reference_sequence, sa, query)
            print(f"Query {i+1}: Positions {positions}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
