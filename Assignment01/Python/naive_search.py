import argparse
import gzip


def read_fasta(file_path):
    """Reads a FASTA file and returns the sequence as a string."""
    sequence = []
    with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):  # Ignore header lines
                sequence.append(line.strip())
    return "".join(sequence)


def naive_search(sequence, query):
    """Performs naive search to find all occurrences of the query in the sequence."""
    positions = []
    start = 0
    while True:
        start = sequence.find(query, start)
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
        reference_sequence = read_fasta(args.reference)
        query_sequences = read_fasta(args.query).split()  # Assuming queries are on separate lines

        # Ensure query count does not exceed available queries
        if args.query_ct > len(query_sequences):
            print(f"Error: The specified query count ({args.query_ct}) exceeds the number of queries in the file.")
            return

        # Perform naive search for the specified number of queries
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
