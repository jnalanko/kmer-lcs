#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

sdsl::int_vector<> build_lcs_superalphabet_algorithm(const sbwt::plain_matrix_sbwt_t& SBWT){
    const sdsl::bit_vector& A_bits = SBWT.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = SBWT.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = SBWT.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = SBWT.get_subset_rank_structure().T_bits;
    int64_t k = SBWT.get_k();
    int64_t n_nodes = SBWT.number_of_subsets();

    return sdsl::int_vector<>(); // TODO
}