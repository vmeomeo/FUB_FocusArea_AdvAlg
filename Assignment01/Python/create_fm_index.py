import iv2py as iv
import argparse


def read_fasta(file_path):
    sequence = []
    for record in iv.fasta.reader(file_path):
        sequence.append(record.seq)  # Extract only the sequence part
    return sequence


# Helper function for reading reference genome
def create_fm_index(reference_file):
    index = iv.fmindex(reference=reference_file, samplingRate=16)
    return index




def main():
    parser = argparse.ArgumentParser(description="FM_index search algorithm for DNA sequences.")
    parser.add_argument("-reference", type=str, required=True,
                        help="Path to the file containing the reference sequence (FASTA format).")


    args = parser.parse_args()

    try:
        # Read reference and query sequences
        reference_sequences = read_fasta(args.reference)
        fm_index = create_fm_index(reference_sequences)
        print("index created")
        fm_index.save("./fm_index_partial.idx")
        print("index saved")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
