#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include "peak_rss.hpp"
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

#include "lcs_naive_algorithm.hpp"
#include "lcs_basic_algorithm.hpp"
#include "lcs_basic_parallel_algorithm.hpp"
#include "lcs_linear_algorithm.hpp"
#include "lcs_superalphabet_algorithm.hpp"

#include "SeqIO.hh"


using namespace std;
using namespace sbwt;


int64_t get_microseconds(){
    // Get the current time in microseconds
    auto start = chrono::high_resolution_clock::now();
    int64_t microseconds = chrono::duration_cast<chrono::microseconds>(start.time_since_epoch()).count();
    return microseconds;
}

// Run an algorithm and print the time it took to run
sdsl::int_vector<> run_with_timing(const sbwt::plain_matrix_sbwt_t& sbwt, std::function<sdsl::int_vector<>(const sbwt::plain_matrix_sbwt_t&)> f, const string& algorithm_name, ostream& out){
    
    cerr << "Running " << algorithm_name << endl;

    int64_t start = get_microseconds();
    sdsl::int_vector<> lcs = f(sbwt);
    int64_t end = get_microseconds();

    out << algorithm_name << "," << (end-start) / 1e6 << "," << getPeakRSS() / (double)(1<<20) << "," << sbwt.get_k() << "," << sbwt.number_of_subsets() << endl;
    return lcs;
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Usage: " << argv[0] << " <sbwt-file> <algorithm> <output-file>" << endl;
        cerr << "The sbwt file must be in plain-matrix format." << endl;
        cerr << "The algorithm is one of: naive, basic, basic-parallel-x superalphabet-x, linear," << endl;
        cerr << "where the x in superalphabet-x is the number of characters in each supercharacter," << endl;
        cerr << "and the x in basic-parallel-x is the number of threads." << endl;
        cerr << "Appends one line to the output file containing the result," << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);
    string algorithm = string(argv[2]);
    string outfile = string(argv[3]);
    

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    sbwt::throwing_ofstream out(outfile, ios::app); // Append mode
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
        return 1;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    sdsl::int_vector<> LCS;

    if(algorithm == "naive"){
        LCS = run_with_timing(sbwt, lcs_naive_algorithm, "Naive", out.stream);
    } else if(algorithm == "basic"){
        LCS = run_with_timing(sbwt, lcs_basic_algorithm, "Basic", out.stream);
    } else if(algorithm.substr(0, 15) == "basic-parallel-"){
        int64_t num_threads = stoll(algorithm.substr(15));
        auto basic_parallel_algorithm = [num_threads](const sbwt::plain_matrix_sbwt_t& sbwt){
            return lcs_basic_parallel_algorithm(sbwt, num_threads);
        };
        LCS = run_with_timing(sbwt, basic_parallel_algorithm, "Basic-Parallel-" + to_string(num_threads), out.stream);
        

    } else if(algorithm == "linear"){
        LCS = run_with_timing(sbwt, lcs_linear_algorithm, "Linear", out.stream);
    } else if(algorithm.substr(0, 14) == "superalphabet-"){
        int64_t superalphabet_size = stoll(algorithm.substr(14));
        auto superalphabet_algorithm = [superalphabet_size](const sbwt::plain_matrix_sbwt_t& sbwt){
            return lcs_superalphabet_algorithm(sbwt, superalphabet_size);
        };
        LCS = run_with_timing(sbwt, superalphabet_algorithm, "Superalphabet-" + to_string(superalphabet_size), out.stream);
    } else {
        cerr << "Error: unknown algorithm " << algorithm << endl;
        return 1;
    }
    if (argc == 5){
        string lcs_outfile = string(argv[4]);
        std::ofstream LCS_out(lcs_outfile + ".sdsl");
        sdsl::serialize(LCS, LCS_out);
    }
}
