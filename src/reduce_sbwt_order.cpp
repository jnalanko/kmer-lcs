#include <iostream>
#include <string>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

#include "lcs_basic_algorithm.hpp"
#include "reduce_sbwt_order.hpp"

using namespace std;
using namespace sbwt;

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);
    string out_file = string(argv[2]);
    int64_t new_k = stoll(argv[3]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    sbwt::throwing_ofstream out(out_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
        return 1;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;

    sdsl::int_vector lcs = lcs_basic_algorithm(sbwt);
    sbwt::plain_matrix_sbwt_t new_sbwt = reduce_sbwt_order(sbwt, lcs, new_k);

    cerr << "Reduced SBWT has " << new_sbwt.number_of_subsets() << " subsets" << endl;

    serialize_string("plain-matrix", out.stream);
    new_sbwt.serialize(out.stream);

}
