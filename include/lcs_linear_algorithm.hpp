#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

sdsl::int_vector<> lcs_linear_algorithm(const sbwt::plain_matrix_sbwt_t& SBWT){
    const int64_t n_nodes = SBWT.number_of_subsets();
    const sdsl::bit_vector& A_bits = SBWT.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = SBWT.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = SBWT.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = SBWT.get_subset_rank_structure().T_bits;

    const sdsl::rank_support_v5<>& A_bits_rs = SBWT.get_subset_rank_structure().A_bits_rs;
    const sdsl::rank_support_v5<>& C_bits_rs = SBWT.get_subset_rank_structure().C_bits_rs;
    const sdsl::rank_support_v5<>& G_bits_rs = SBWT.get_subset_rank_structure().G_bits_rs;
    const sdsl::rank_support_v5<>& T_bits_rs = SBWT.get_subset_rank_structure().T_bits_rs;

    const sdsl::bit_vector* DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};
    const sdsl::rank_support_v5<>* DNA_rs[4] = {&A_bits_rs, &C_bits_rs, &G_bits_rs, &T_bits_rs};

    const int64_t k = SBWT.get_k();
    const vector<int64_t>& C = SBWT.get_C_array();
    const int sigma = C.size();

    sdsl::int_vector<> lcs(n_nodes, k, 64 - __builtin_clzll(k)); // Enough bits per element to store values from 0 to k
    lcs[0]=lcs[1]=0;
    vector<pair<uint64_t, uint64_t>> I;
    I.push_back({0,n_nodes});

    vector<pair<uint64_t, uint64_t>> _I = {{0,0}}; // $ interval
    for(int64_t i = 0; i < k; i++) {
        cerr << "Round " << i << "/" << k-1 << ", intervals: " << I.size() << endl;
        while (!I.empty()){
            pair<uint64_t, uint64_t> l_r = I.back();
            I.pop_back();
            // Enumerate Right
            for (int c = 0; c < sigma; c++){
                const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[c]);
                const sdsl::rank_support_v5<>& Bit_rs = *(DNA_rs[c]);
                // Extend right
                // if rank < 0 no possible extension
                int64_t l = C[c] + Bit_rs.rank(l_r.first);
                int64_t r = C[c] + Bit_rs.rank(l_r.second + 1) -1;
                int64_t rank = r - l ;
                if (rank >= 0 && r < n_nodes-1 && lcs[r+1] == k){
                    lcs[r+1] = i;
                    _I.push_back({ l, r });
                }
            }
        }
        I = std::move(_I);
        _I.clear();
    }
    return lcs;
}
