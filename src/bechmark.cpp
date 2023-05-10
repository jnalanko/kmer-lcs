#include <iostream>
#include <string>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

#include "lcs_naive_algorithm.hpp"
#include "lcs_basic_algorithm.hpp"
#include "lcs_linear_algorithm.hpp"
#include "lcs_superalphabet_algorithm.hpp"

using namespace std;
using namespace sbwt;

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file" << endl;
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

    cerr << "Building LCS with naive algorithm" << endl;
    sdsl::int_vector naive = lcs_naive_algorithm(sbwt);
    cerr << "Building LCS with basic algorithm" << endl;
    sdsl::int_vector basic = lcs_basic_algorithm(sbwt);
    cerr << "Building LCS with superalphabet algorithm" << endl;
    sdsl::int_vector superalphabet = lcs_superalphabet_algorithm(sbwt,2);
    cerr << "Building LCS with linear algorithm" << endl;
    sdsl::int_vector linear = lcs_linear_algorithm(sbwt);

    if(basic != superalphabet) cerr << "wrong answer: basic != superalphabet" << endl;
    if(superalphabet != linear) cerr << "wrong answer: superalphabet != linear" << endl;

}
