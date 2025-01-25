import argparse
import gzip
#Willll it work?
def read_fasta(file_path):
    """Reads sequences from a FASTA file (supports gzip)."""
    sequences = []
    with gzip.open(file_path, 'rt') as f:
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if sequence:  # If there's a sequence collected, add it
                    sequences.append("".join(sequence))
                    sequence = []  # Reset for the next sequence
            else:
                sequence.append(line)  # Append sequence line
        if sequence:  # Add the last sequence if present
            sequences.append("".join(sequence))
    return sequences


def naive_search(reference, query):
    """Performs naive search to find all occurrences of the query in the reference."""
    positions = []
    start = 0
    while True:
        start = reference.find(query, start)
        if start == -1:  # No more matches
            break
        positions.append(start)
        start += 1
    return positions

def main():
    parser = argparse.ArgumentParser(description="Naive search algorithm for DNA sequences.")
    parser.add_argument("-reference", type=str, required=True,
                        help="Path to the file containing the reference sequence (FASTA format).")
    parser.add_argument("-query", type=str, required=True,
                        help="Path to the file containing the query sequences (FASTA format).")
    parser.add_argument("-query_ct", type=int, required=True,
                        help="Number of query sequences to search for.")

    args = parser.parse_args()

    try:
        # Read reference and query sequences
        reference_sequences = read_fasta(args.reference)
        query_sequences = read_fasta(args.query)

        # Combine all reference sequences into a single string
        reference_sequence = "".join(reference_sequences)

        # Ensure query count does not exceed available queries
        if args.query_ct > len(query_sequences):
            print(f"Error: The specified query count ({args.query_ct}) exceeds the number of queries available ({len(query_sequences)}).")
            return

        # Perform naive search for each query sequence
        for i in range(args.query_ct):
            query = query_sequences[i]
            positions = naive_search(reference_sequence, query)
            print(f"Positions: {positions}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
