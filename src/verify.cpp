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
#include "lcs_basic_parallel_algorithm.hpp"

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
        cerr << "Usage: " << argv[0] << " <sbwt-file>" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
        return 1;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    sdsl::int_vector naive = lcs_naive_algorithm(sbwt);
    sdsl::int_vector basic = lcs_basic_algorithm(sbwt);
    sdsl::int_vector basic_parallel = lcs_basic_parallel_algorithm(sbwt);
    sdsl::int_vector superalphabet_2 = lcs_superalphabet_algorithm(sbwt, 2);
    sdsl::int_vector superalphabet_4 = lcs_superalphabet_algorithm(sbwt, 4);
    sdsl::int_vector linear = lcs_linear_algorithm(sbwt);

    //if(basic != naive) cerr << "Basic and naive algorithms do not agree" << endl;
    if(basic != basic_parallel) cerr << "Basic and basic parallel algorithms do not agree" << endl;
    //if(basic != superalphabet_2) cerr << "Basic and superalphabet-2 algorithms do not agree" << endl;
    //if(basic != superalphabet_4) cerr << "Basic and superalphabet-4 algorithms do not agree" << endl;
    if(basic != linear) cerr << "Basic and linear algorithms do not agree" << endl;
    //else cerr << "All algorithms agree" << endl;

}