#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

#include "lcs_naive_algorithm.hpp"
#include "lcs_basic_algorithm.hpp"
#include "lcs_linear_algorithm.hpp"
#include "lcs_superalphabet_algorithm.hpp"

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

    out << algorithm_name << "," << (end-start) / 1e6 << endl;
    return lcs;
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file and the name of an output file" << endl;
        cerr << "Usage: " << argv[0] << " <sbwt-file> <output-file>" << endl;
        cerr << "The output file will be a CSV file with two columns: algorithm name and running time in seconds" << endl;
        cerr << "The output file will be overwritten if it already exists" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);
    string outfile = string(argv[2]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    sbwt::throwing_ofstream out(outfile);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
        return 1;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    auto superalphabet_algorithm_2 = [](const sbwt::plain_matrix_sbwt_t& sbwt){
        return lcs_superalphabet_algorithm(sbwt, 2);
    };

    auto superalphabet_algorithm_4 = [](const sbwt::plain_matrix_sbwt_t& sbwt){
        return lcs_superalphabet_algorithm(sbwt, 4);
    };

    sdsl::int_vector naive = run_with_timing(sbwt, lcs_naive_algorithm, "Naive", out.stream);
    sdsl::int_vector basic = run_with_timing(sbwt, lcs_basic_algorithm, "Basic", out.stream);
    sdsl::int_vector superalphabet_2 = run_with_timing(sbwt, superalphabet_algorithm_2, "Superalphabet-2", out.stream);
    sdsl::int_vector superalphabet_4 = run_with_timing(sbwt, superalphabet_algorithm_4, "Superalphabet-4", out.stream);
    sdsl::int_vector linear = run_with_timing(sbwt, lcs_linear_algorithm, "Linear", out.stream);

    if(naive != basic) cerr << "Naive and basic algorithms do not agree" << endl;
    if(naive != superalphabet_2) cerr << "Naive and superalphabet-2 algorithms do not agree" << endl;
    if(naive != superalphabet_4) cerr << "Naive and superalphabet-4 algorithms do not agree" << endl;
    if(naive != linear) cerr << "Naive and linear algorithms do not agree" << endl;

}
