#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

using namespace std;

// This is a supporting function for reduce_sbwt_order
sdsl::bit_vector to_sdsl(const vector<bool>& vec){
    sdsl::bit_vector bv(vec.size(), 0);
    for(int64_t i = 0; i < vec.size(); i++) bv[i] = vec[i];
    return bv;
}

// This is a supporting function for reduce_sbwt_order
// Merges k-mers that have longest common suffix length of at least new_k
// The resulting SBWT may have redundant dummies, and the number_of_kmers variable is not set.
sbwt::plain_matrix_sbwt_t merge_kmers(const sbwt::plain_matrix_sbwt_t& SBWT, const sdsl::int_vector<>& LCS, int64_t new_k){
    vector<const sdsl::bit_vector*> matrix = {
        &SBWT.get_subset_rank_structure().A_bits,
        &SBWT.get_subset_rank_structure().C_bits,
        &SBWT.get_subset_rank_structure().G_bits,
        &SBWT.get_subset_rank_structure().T_bits};

    int64_t n_nodes = SBWT.number_of_subsets();
    int64_t sigma = matrix.size();
    vector<char> alphabet;
    for(char c = 0; c < sigma; c++) alphabet.push_back(c);

    vector<vector<bool>> new_matrix(4);
    vector<bool> new_suffix_group_starts;

    // Merge old k-mers having the same suffix of length equal to the new k
    for(int64_t i = 0; i < n_nodes; i++){
        if(LCS[i] >= new_k){
            // Add bits to the previous column
            for(char c : alphabet)
                new_matrix[c].back() = new_matrix[c].back() || (*matrix[c])[i];
        } else{
            
            new_suffix_group_starts.push_back(LCS[i] < new_k-1);
            
            // Copy the old column
            for(char c : alphabet) new_matrix[c].push_back((*matrix[c])[i]);
        }
    }

    // Zero out columns that are not at the starts of suffix groups
    int64_t prev_suffix_group_start = 0;
    for(int64_t i = 0; i < new_suffix_group_starts.size(); i++){
        if(new_suffix_group_starts[i] == 0){
            for(char c : alphabet){
                new_matrix[c][prev_suffix_group_start] = new_matrix[c][prev_suffix_group_start] | new_matrix[c][i];
                new_matrix[c][i] = 0;
            }
        } else prev_suffix_group_start = i;
    }

    // This SBWT is functional, but can still have redundant dummies
    return sbwt::plain_matrix_sbwt_t(
        to_sdsl(new_matrix[0]),
        to_sdsl(new_matrix[1]),
        to_sdsl(new_matrix[2]),
        to_sdsl(new_matrix[3]),
        to_sdsl(new_suffix_group_starts),
        new_k, 0, 0); // number of k-mers is set to 0 because we don't know it, but that does not matter

}

// This is a supporting function for reduce_sbwt_order
void filter_bit_vectors(vector<sdsl::bit_vector*> vecs, const sdsl::bit_vector* to_delete){
    int64_t j = 0;
    for(int64_t i = 0; i < to_delete->size(); i++){
        if((*to_delete)[i] == false){ // keep
            for(sdsl::bit_vector* vec : vecs){
                (*vec)[j] = (*vec)[i];
            }
            j++;
        }
    }
    for(sdsl::bit_vector* vec : vecs){
        vec->resize(j);
    }
}


// This is a supporting function for reduce_sbwt_order
// Returns modified four sbwt bit vectors and the new streaming support
std::array<sdsl::bit_vector, 5> delete_dummies(const sbwt::plain_matrix_sbwt_t& SBWT){
    int64_t old_n = SBWT.number_of_subsets();
    sdsl::bit_vector delete_marks(old_n, 0);
    sdsl::bit_vector A = SBWT.get_subset_rank_structure().A_bits;
    sdsl::bit_vector C = SBWT.get_subset_rank_structure().C_bits;
    sdsl::bit_vector G = SBWT.get_subset_rank_structure().G_bits;
    sdsl::bit_vector T = SBWT.get_subset_rank_structure().T_bits;
    sdsl::bit_vector ss = SBWT.get_streaming_support();

    string ACGT = "ACGT";
    std::function<bool(int64_t, int64_t)> dfs = [&](int64_t v, int64_t depth){
        if(depth == SBWT.get_k()-1){ // Recursion base case
            // Return true iff the suffix group of v has size greater than 1
            bool to_delete = v < old_n-1 && SBWT.get_streaming_support()[v+1] == 0;
            if(to_delete){
                // Push the outgoing edges to the next member of the suffix group
                A[v+1] = A[v]; A[v] = 0;
                C[v+1] = C[v]; C[v] = 0;
                G[v+1] = G[v]; G[v] = 0;
                T[v+1] = T[v]; T[v] = 0;
                ss[v+1] = ss[v]; ss[v] = 0; // The new start of the suffix group
            }
            return to_delete;
        } else{ // Push children
            bool all_children_deleted = true;
            for(char c : ACGT){
                int64_t u = SBWT.forward(v, c);
                if(u != -1){
                    bool to_delete = dfs(u, depth+1);
                    if(to_delete){
                        delete_marks[u] = 1; // Delete node

                        // Delete edge to the node
                        if (c == 'A') A[v] = 0;
                        else if(c == 'C') C[v] = 0;
                        else if(c == 'G') G[v] = 0;
                        else if(c == 'T') T[v] = 0;
                    }
                    all_children_deleted &= to_delete;
                }
            }            
            return all_children_deleted;
        }
    };

    dfs(0,0);
    filter_bit_vectors({&A, &C, &G, &T, &ss}, &delete_marks);
    return {A,C,G,T,ss};
}

// Return a SBWT with the new k
sbwt::plain_matrix_sbwt_t reduce_sbwt_order(const sbwt::plain_matrix_sbwt_t& SBWT, const sdsl::int_vector<>& LCS, int64_t new_k){

    cerr << "Merging k-mers" << endl;
    sbwt::plain_matrix_sbwt_t merged_SBWT = merge_kmers(SBWT, LCS, new_k);
    // This SBWT is functional, but can still have redundant dummies

    cerr << "Merged SBWT has " << merged_SBWT.number_of_subsets() << " subsets" << endl;

    cerr << "Deleting redundant dummies" << endl;
    auto [final_A, final_C, final_G, final_T, final_suffix_groups] = delete_dummies(merged_SBWT);

    cerr << "Building the final SBWT" << endl;
    return sbwt::plain_matrix_sbwt_t(
        final_A, final_C, final_G, final_T, 
        final_suffix_groups,
        new_k, 0, 0
    );
}
