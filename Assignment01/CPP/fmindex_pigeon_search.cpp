#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream) {
        reference.push_back(record.sequence());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};
    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text
    // Pseudo code (might be wrong):
    // for query in queries:
    //      parts[3] = cut_query(3, query);
    //      for p in {0, 1, 2}:
    //          for (pos in search(index, part[p]):
    //              if (verify(ref, query, pos +- ....something)):
    //                  std::cout << "found something\n"
    for (auto const& query : queries) {
        size_t query_length = query.size();
        size_t part_length = query_length / 3;

        // Split the query into 3 parts
        //auto part1 = query | seqan3::views::slice(0, part_length);
        //auto part2 = query | seqan3::views::slice(part_length, 2 * part_length);
        //auto part3 = query | seqan3::views::slice(2 * part_length, query_length);
	//
	//
        std::vector<seqan3::dna5> part1(query.begin(), query.begin() + query.size() / 3);
        std::vector<seqan3::dna5> part2(query.begin() + query.size() / 3, query.begin() + 2 * query.size() / 3);
        std::vector<seqan3::dna5> part3(query.begin() + 2 * query.size() / 3, query.end());

        std::vector<std::vector<seqan3::dna5>> parts{part1, part2, part3};

        for (size_t part_idx = 0; part_idx < 3; ++part_idx) {
            auto search_results = seqan3::search(parts[part_idx], index, cfg);

            for (auto const& result : search_results) {
                size_t pos = result.reference_begin_position();

                // Verify the full query match
                size_t start = (part_idx == 0) ? pos : pos - part_idx * part_length;
                if (start + query_length <= reference[0].size()) { // Ensure valid range
                    auto ref_subseq = reference[0] | seqan3::views::slice(start, start + query_length);

                    if (std::equal(query.begin(), query.end(), ref_subseq.begin())) {
                        seqan3::debug_stream << "Query found at position: " << start << "\n";
                    }
                }
            }
        }
    }
    return 0;
}
