import iv2py as iv


def read_fasta(file_path):
    sequence = []
    for record in iv.fasta.reader(file_path):
        sequence.append(record.seq)  # Extract only the sequence part
    return sequence
reference= read_fasta("../data/hg38_partial.fasta.gz")
#print(len(reference[0])) 100M

# Helper function for reading reference genome
def create_fm_index(reference_file):
    index = iv.fmindex(reference=reference_file, samplingRate=16)
    return index
fm_in= create_fm_index(reference)
#fm_index=create_fm_index(reference)

#
# # Helper function for searching queries
# def search_queries(index, query_file, output_file):
#     """
#     Search for query sequences in the FM-index and save results to an output file.
#     """
#     with open(output_file, "w") as out_file:
#         for record in iv.fasta.reader(file=query_file):
#             res = index.search(record.seq)
#             out_file.write(f"Marker '{record.seq}' found at positions: {res}\n")
#     print(f"Results saved to {output_file}")
#
# # Helper function for benchmarking
# def benchmark_function(function, *args, **kwargs):
#     """
#     Benchmark runtime and memory usage of a function.
#     """
#     start_time = time.time()
#     memory_before = memory_usage()[0]
#     function(*args, **kwargs)
#     memory_after = memory_usage()[0]
#     print(f"Runtime: {time.time() - start_time:.2f} seconds")
#     print(f"Memory used: {memory_after - memory_before:.2f} MB")
#
# # Main execution
# def main(reference_file, query_file, fmindex_file="fmindex", output_file="search_results.txt"):
#     # Load reference genome and build/load FM-index
#     index = load_reference(reference_file, fmindex_file)
#
#     # Benchmark search process
#     print("Starting query search...")
#     benchmark_function(search_queries, index, query_file, output_file)
#
# # Specify file paths
# reference_file = "hg38_partial.fasta.gz"
# query_file = "illumina_reads_40.fasta.gz"
# fmindex_file = "fmindex"
# output_file = "search_results.txt"
#
# # Run the main program
# main(reference_file, query_file, fmindex_file, output_file)