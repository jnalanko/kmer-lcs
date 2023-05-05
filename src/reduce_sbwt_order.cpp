#include <iostream>
#include <string>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

#include "lcs_basic_algorithm.hpp"
#include "lcs_linear_algorithm.hpp"
#include "reduce_sbwt_order.hpp"

using namespace std;
using namespace sbwt;

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Usage: " << argv[0] << " index.sbwt output_prefix 16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 255" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);
    string out_prefix = string(argv[2]);
    vector<int64_t> new_k_values;
    for(int64_t i = 3; i < argc; i++){
        int64_t k = stoll(argv[i]);
        new_k_values.push_back(k);
    }

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
        return 1;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    int64_t max_k = *std::max_element(new_k_values.begin(), new_k_values.end());
    if(max_k > sbwt.get_k()){
        cerr << "Error: new k must be smaller or equal to the original. New = " << max_k << ", original = " << sbwt.get_k() << endl;
        return 1;
    }

    cerr << "Building the LCS" << endl;
    sdsl::int_vector lcs = lcs_linear_algorithm(sbwt);

    for(int64_t k : new_k_values){
        cerr << "Reducing to order " << k << endl;
        sbwt::throwing_ofstream out_file(out_prefix + "-" + to_string(k) + ".sbwt", ios::binary);
        
        sbwt::plain_matrix_sbwt_t new_sbwt = reduce_sbwt_order(sbwt, lcs, k);

        cerr << "SBWT of order " << k << " has " << new_sbwt.number_of_subsets() << " subsets" << endl;

        serialize_string("plain-matrix", out_file.stream);
        new_sbwt.serialize(out_file.stream);
    }
}
