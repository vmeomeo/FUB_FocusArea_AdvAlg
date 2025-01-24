import iv2py as iv
import argparse
import pickle

def read_fasta(file_path):
    sequence = []
    for record in iv.fasta.reader(file_path):
        sequence.append(record.seq)  # Extract only the sequence part
    return sequence


def find_match(query,index):
    result = index.search(query)
    return result



def load_fm_index(index_path):
    """Load the FM-index from a file."""
    with open(index_path, "rb") as f:
        return pickle.load(f)

def main():
    parser = argparse.ArgumentParser(description="FM_index search algorithm for DNA sequences.")
    parser.add_argument("-index", type=str, required=True,
                        help="Path to the file containing the fm_index.")
    parser.add_argument("-query", type=str, required=True,
                        help="Path to the file containing the query sequences (FASTA format).")
    parser.add_argument("-query_ct", type=int, required=True,
                        help="Number of query sequences to search for.")

    args = parser.parse_args()

    try:
        # Read reference and query sequences
        query_sequences = read_fasta(args.query)



        # Ensure query count does not exceed available queries
        if args.query_ct > len(query_sequences):
            print(f"Error: The specified query count ({args.query_ct}) exceeds the number of queries available ({len(query_sequences)}).")
            return
        fm_index=load_fm_index(args.index)
        # Perform naive search for each query sequence
        for i in range(args.query_ct):
            query = query_sequences[i]
            positions= find_match(query,fm_index)

            print(f"Positions: {positions}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()