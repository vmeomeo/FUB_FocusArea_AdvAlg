import iv2py as iv
import argparse

def read_fasta(file_path):
    sequence = []
    for record in iv.fasta.reader(file_path):
        sequence.append(record.seq)  # Extract only the sequence part
    return sequence



def find_match(query,index,errors):
    result = index.search(query, k=errors)
    return result



def load_fm_index(index_path):
    """Load the FM-index from a file."""
    index = iv.fmindex(path=index_path)
    return index

def main():
    parser = argparse.ArgumentParser(description="FM_index search algorithm for DNA sequences.")
    parser.add_argument("-index", type=str, required=True,
                        help="Path to the file containing the fm_index.")
    parser.add_argument("-query", type=str, required=True,
                        help="Path to the file containing the query sequences (FASTA format).")
    parser.add_argument("-errors", type=int, required=True,
                        help="Number of errors in the sequences.")

    args = parser.parse_args()

    try:
        # Read reference and query sequences
        query_sequences = read_fasta(args.query)
        print("Query read")
        fm_index=load_fm_index(args.index)
        print("Index loaded")
        # Perform naive search for each query sequence
        for i in range(10000):
            query = query_sequences[i]
            positions= find_match(query,fm_index,args.errors)

            print(f"Positions: {positions}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
