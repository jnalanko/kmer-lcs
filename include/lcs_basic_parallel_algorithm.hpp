#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include <vector>

using namespace std;

// End is one past end
void lcs_worker_thread(char** output_starts, const char* input, const sdsl::bit_vector** bits, int64_t start, int64_t end){
    for(int64_t i = start; i < end; i++){
        for(char sigma = 0; sigma < 4; sigma++){
            if((*(bits[sigma]))[i] == 1){
                *(output_starts[sigma]) = input[i];
                output_starts[sigma]++;
            }
        }
    }
}

sdsl::int_vector<> lcs_basic_parallel_algorithm(const sbwt::plain_matrix_sbwt_t& SBWT){

    const sdsl::bit_vector& A_bits = SBWT.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = SBWT.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = SBWT.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = SBWT.get_subset_rank_structure().T_bits;
    int64_t k = SBWT.get_k();
    int64_t n_nodes = SBWT.number_of_subsets();

    vector<int64_t> C_array(4);

    vector<char> last; // last[i] = incoming character to node i
    last.push_back('$');

    C_array[0] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(A_bits[i]) last.push_back('A');

    C_array[1] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(C_bits[i]) last.push_back('C');

    C_array[2] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(G_bits[i]) last.push_back('G');
    
    C_array[3] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(T_bits[i]) last.push_back('T');

    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }

    sdsl::bit_vector mismatch_found_marks(n_nodes, 0);
    sdsl::int_vector<> lcs(n_nodes, 0, 64 - __builtin_clzll(k-1)); // Enough bits per element to store values from 0 to k-1

    for(int64_t round = 0; round < k; round++){
        cerr << "Round " << round << "/" << k-1 << endl;
        for(int64_t i = 0; i < n_nodes; i++){
            if(mismatch_found_marks[i] == 0 && (i == 0 || last[i] != last[i-1])){
                mismatch_found_marks[i] = 1;
                lcs[i] = round;
            }
        }

        // Propagate the labels one step forward in the graph
        vector<char> propagated(n_nodes, '$');
        char* A_ptr = propagated.data() + C_array[0];
        char* C_ptr = propagated.data() + C_array[1];
        char* G_ptr = propagated.data() + C_array[2];
        char* T_ptr = propagated.data() + C_array[3];
        char* ptrs[4] = {A_ptr, C_ptr, G_ptr, T_ptr};
        const sdsl::bit_vector* bits[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
        lcs_worker_thread(ptrs, last.data(), bits, 0, n_nodes);
        last = propagated;
    }

    return lcs;
}