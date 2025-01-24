import iv2py as iv
import os
import gzip
from Bio import SeqIO
import time
import argparse

# Load reference genome
def load_reference(file_path):
    """ Load a fasta.gz file
        return combined sequence as a string """  
    with gzip.open(file_path, 'rt') as f:
        sequences = [str(record.seq) for record in SeqIO.parse(f, 'fasta')]
    return ''.join(sequences) #if more than one sequence is given, combine them into one string

# Load marker files
def load_markers(file_path):
    """ Load marker sequences from fasta.gz file
        return marker sequences as individual strings """
    with gzip.open(file_path, 'rt') as f:
        return [str(record.seq) for record in SeqIO.parse(f, 'fasta')]
    
# Search the markers in the reference genome
def search_suffix_array(suffix_array, reference, marker):
    """ Find all occurences of a marker in reference genome using mlr"""
    left, right = 0, len(suffix_array) - 1
    start, end = -1, -1

    # Left boundary search
    while left <= right:
        mid = (left + right) // 2
        suffix = reference[suffix_array[mid]:suffix_array[mid] + len(marker)]
        if suffix < marker:
            left = mid + 1
        elif suffix > marker:
            right = mid - 1
        else:
            start = mid
            right = mid - 1 # Continue search for the start of the range

    # No match found
    if start == -1:
        return []
    
    # Right boundary search
    left, right = start, len(suffix_array) - 1
    while left <= right:
        mid = (left + right) // 2
        suffix = reference[suffix_array[mid]:suffix_array[mid] + len(marker)]
        if suffix == marker: 
            end = mid 
            left = mid + 1 # Countinue search for the end of the range
        else:
            right = mid - 1

    # Collect and return all positions within the range
    return suffix_array[start:end + 1]


    




def main(path_reference_genome, path_query_file, num_query):
    start_time = time.time()
    
    reference_genome = load_reference(path_reference_genome)
    print('Reference genome loaded.')

    suffix_array = iv.create_suffixarray(reference_genome)
    print("Suffix array built.")


    markers = load_markers(path_query_file)
    print('Queries loaded.')

    # Check if the number of queries to search for exceeds the number of queries available
    if num_query > len(markers):
        markers = markers * (num_query // len(markers)) + markers[:num_query % len(markers)]
    else:
        markers = markers[:num_query]

    # Find all markers in the reference genome
    results = {}
    for marker in markers:
        position = search_suffix_array(suffix_array, reference_genome, marker)
        results[marker] = position

    search_time = time.time() - start_time
    

    # Print the results
    for marker, position in results.items():
        if position:
            print(f"Marker '{marker}' found at position {position}.")

    print(f"Total search time: {search_time:.2f} seconds.")

if __name__ == "__main__":
    # Set up argparse to handle command line arguments
    parser = argparse.ArgumentParser(description="Search for markers in a reference genome using suffix array search.")
    
    # Define arguments
    parser.add_argument('-r', '--reference', type=str, required=True, help="Path to the reference genome file (fasta.gz).")
    parser.add_argument('-q', '--query', type=str, required=True, help="Path to the query markers file (fasta.gz).")
    parser.add_argument('-q_c', '--query_count', type=int, required=True, help="Number of queries to search for.")
    
    # Parse arguments
    args = parser.parse_args()

    # Call the main function with arguments
    main(args.reference, args.query, args.query_count)