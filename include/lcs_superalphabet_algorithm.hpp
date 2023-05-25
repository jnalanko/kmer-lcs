#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include <cmath>
#include "bitset"

int64_t get_char_idx(char c){
    switch(c){
        case '$': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        default: return 0;
    }
}

char get_char_from_key(char k){
    if (k == 0){ return 0;}
    char c1 = (k - 1) % 5;// 0=$,1=A,2=C,3=G,4=T
    char c2 = ceil(static_cast<double>(k)/5) ;
    char super_char =  (c1 << 3) | c2;
    return super_char;
}

char get_key_from_char(char c1, char c2){
    int k = c2 * 5 - (4 - c1);
    return k;
}

char get_keyab_from_key(char key, const vector<char>& alphabet){
    auto it = std::find(alphabet.begin(), alphabet.end(), key);
    if (it != alphabet.end()) {
        return it - alphabet.begin();
    }
    return  - 1;
}

char get_keyab_from_char(char c1, char c2, const vector<char>& alphabet){
    int k = c2 * 5 - (4 - c1);
    return get_keyab_from_key(k, alphabet);
}

char get_key_from_keyab(char keyab, const vector<char>& alphabet){
    char key = alphabet[keyab];
    return key;
}

struct concat_str{
    vector<char> concat_v;
    sdsl::bit_vector concat_b;
};

concat_str get_super_concat(char base, const vector<char>& concat_vector, const sdsl::bit_vector& concat_bits, const vector<int64_t>& C_array, const vector<char>& super_alphabet, const vector<char>& alphabet){
    concat_str result;
    vector<char> res_concat_v;
    vector<bool> bv;
    sdsl::bit_vector::select_1_type select_concat_bits(&concat_bits);

    // Todo remove select
    // Scan concat_bv once and add a pointer for every element in the alphabet
    vector<int64_t> select_v = {0};
    int64_t i=1, subset_counter = 0; // In the first pos there is a 1 = {} = $$...
    for (int64_t k = 0; k < C_array.size(); k++){
        while (i < concat_bits.size()){
            if (concat_bits[i]){ subset_counter++;}
            if(subset_counter == C_array[k] ){
                select_v.push_back(i);
                i++;
                break;
            }
            i++;
        }
    }

    vector<int64_t> C_copy = {0};
    C_copy.insert(C_copy.end(), C_array.begin(), C_array.end());
    int64_t subsets = 0;
    for (int64_t x = 0; x < concat_vector.size(); x++){
        if (concat_bits[x+subsets]){
            subsets++;
            bv.push_back(1);
            // (1) Check you are not starting from an empty subset -> empty
            if (concat_bits[x+subsets]){
                x--;
                continue; }
        }
        // Find the start of the curr_char-block and append to curr_char the values found in the subset at the start of its block
        char curr_char = concat_vector[x]; // C1 => curr_char + c2
        C_copy[curr_char]++;
        int64_t x_subset = C_copy[curr_char];

        // Start of the subset of curr_char
        // INDEX IN BITVECTOR
        int64_t x_subset_start_bv = select_v[curr_char] + 1; // later +1 as we want zeros
        // INDEX CONCAT VECTOR
        int64_t x_subset_start_v = x_subset_start_bv - x_subset;//rank 0
        // Scan the subset and add each character found to curr_char
        // (2) Check you are not appending an empty set
        while (x_subset_start_bv < concat_bits.size() && !concat_bits[x_subset_start_bv]){
            char key_ab = get_keyab_from_char(curr_char, concat_vector[x_subset_start_v], super_alphabet);
            res_concat_v.push_back(key_ab);
            x_subset_start_bv ++;
            x_subset_start_v ++;
        }
        select_v[curr_char]=x_subset_start_bv;

        // add necessary 0s to the bitvector
        while(bv.size() < (res_concat_v.size() + subsets)) {bv.push_back(0);}
    }
    sdsl::bit_vector res_concat_b(bv.size());
    for (size_t i = 0; i < bv.size(); i++) {
        res_concat_b[i] = bv[i];
    }
    result.concat_v = res_concat_v;
    result.concat_b = res_concat_b;
    return result;
}


sdsl::int_vector<> lcs_superalphabet_algorithm(const sbwt::plain_matrix_sbwt_t& SBWT, int64_t chars_per_super_alpha_char){
    const sdsl::bit_vector& A_bits = SBWT.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = SBWT.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = SBWT.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = SBWT.get_subset_rank_structure().T_bits;
    int64_t k = SBWT.get_k();
    int64_t n_nodes = SBWT.number_of_subsets();
    vector<int64_t> C_array(4);

    vector<char> last = {0};
    C_array[0] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(A_bits[i]) last.push_back(1);
    C_array[1] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(C_bits[i]) last.push_back(2);
    C_array[2] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(G_bits[i]) last.push_back(3);
    C_array[3] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(T_bits[i]) last.push_back(4);
    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }
    sdsl::bit_vector mismatch_found_marks(n_nodes, 0);
    mismatch_found_marks[0]=1;
    sdsl::int_vector<> lcs(n_nodes, 0, 64 - __builtin_clzll(k-1)); // Enough bits per element to store values from 0 to k-1

    // Propagate the labels one step forward in the graph
    vector<char> propagated(n_nodes, 0);
    int64_t A_ptr = C_array[0];
    int64_t C_ptr = C_array[1];
    int64_t G_ptr = C_array[2];
    int64_t T_ptr = C_array[3];
    for(int64_t i = 0; i < n_nodes-1; i++){
        if(mismatch_found_marks[i+1] == 0 &&  last[i] != last[i+1]){ // last[i-1] has been modified already
            mismatch_found_marks[i+1] = 1;
            //lcs[i+1] = 0; already 0
        }
        if(A_bits[i]) propagated[A_ptr++] = last[i];
        if(C_bits[i]) propagated[C_ptr++] = last[i];
        if(G_bits[i]) propagated[G_ptr++] = last[i];
        if(T_bits[i]) propagated[T_ptr++] = last[i];
    }
    if(A_bits[n_nodes-1]) propagated[A_ptr++] = last[n_nodes-1];
    if(C_bits[n_nodes-1]) propagated[C_ptr++] = last[n_nodes-1];
    if(G_bits[n_nodes-1]) propagated[G_ptr++] = last[n_nodes-1];
    if(T_bits[n_nodes-1]) propagated[T_ptr++] = last[n_nodes-1];

    // CONCAT with the basic alphabet
    // we already have the C array
    vector<char> concat_v;
    vector<bool> bv;
    for(int64_t i = 0; i < n_nodes; i++){
        if(A_bits[i]) { concat_v.push_back(1);}
        if(C_bits[i]) { concat_v.push_back(2);}
        if(G_bits[i]) { concat_v.push_back(3);}
        if(T_bits[i]) { concat_v.push_back(4);}
        bv.push_back(1);
        while (bv.size() < (concat_v.size() + i + 1)){
            bv.push_back(0);
        }
    }
    sdsl::bit_vector concat_b(bv.size());
    for (size_t i = 0; i < bv.size(); i++) {
        concat_b[i] = bv[i];
    }

    vector<char> alphabet = {0,1,2,3,4};//$,A,C,G,T // TODO better alphabet mapping

    // create a super alphabet and 2-C_array
    last = {0};
    vector<char> super_alphabet = {0};
    vector<int64_t> super_C_array = {0};
    for ( size_t i = 1; i < n_nodes; i++ ) { // start from 1 as later i-1 and $$ already ok
        // mismatch in the 1st char
        char key = get_key_from_char(propagated[i], last[i]);
        if (mismatch_found_marks[i]){
            super_alphabet.push_back(key);
            super_C_array.push_back(i);
            // mismatch in the 2nd char
        } else if (!mismatch_found_marks[i] && propagated[i] != propagated[i-1]){
            mismatch_found_marks[i] = 1;
            lcs[i] = 1;
            super_C_array.push_back(i);
            super_alphabet.push_back(key);
            // no mismatch yet
        }
        // Create the first 2 columns according to the new super_C_array
        last.push_back(get_keyab_from_key(key, super_alphabet));

    }
    super_C_array.erase(super_C_array.begin());// remove the first 0

    // CONCAT 2-super alphabet
    auto[super_concat_v, super_concat_b] = get_super_concat(2,concat_v, concat_b, C_array, super_alphabet, alphabet );

// -------------------------------------------------------------------
    // Dump k-mers with 2-super_alphabet
    // Propagate the labels one step forward in the graph

    for(int64_t round = chars_per_super_alpha_char; round < k; round+=chars_per_super_alpha_char){ // works w chars_per_super_alpha_char = 2
        // Propagate the labels one step forward in the graph
        vector<int64_t> super_C_copy = {0};
        super_C_copy.insert(super_C_copy.end(), super_C_array.begin(), super_C_array.end());
        int64_t subsets = 0;
        for(int64_t i = 0; i < super_concat_v.size(); i++){
            if (super_concat_b[i+subsets]){
                subsets++;
                if (super_concat_b[i+subsets]){ // empty subset
                    i--;
                    continue;
                }
            }
            char c = super_concat_v[i];
            propagated[super_C_copy[c]++] = last[subsets-1];
        }

        for(int64_t j = 1; j < n_nodes; j++){
            if(!mismatch_found_marks[j]) {
                char pi = get_char_from_key(get_key_from_keyab(propagated[j], super_alphabet));
                char pi_ = get_char_from_key(get_key_from_keyab(propagated[j-1], super_alphabet));
                // Compare pi and pi-1 with xor and mask the result to see only the last 3 bits
                char compare_res = pi ^ pi_;
                if (compare_res != 0){
                    mismatch_found_marks[j] = 1;
                    // Check only c2 (last 3 bits)
                    compare_res = compare_res & 0x07;
                    if ( compare_res == 0){ // No mismatch in c2 -> mismatch in c1
                        lcs[j] = round + 1;
                    } else { // Mismatch in c2
                        lcs[j] = round;
                    }
                }
            }
        }
        last = propagated;
    }

//    // CONCAT with the 4-super alphabet
//    if(chars_per_super_alpha_char == 4 ){
//        super_ab.clear();
//        //TODO create a super alphabet of 4
//        //todo define empty  char
//        auto[super_concat_v, super_concat_b] = get_super_concat(super_concat_v, super_concat_b, super_alphabet);
//    }

    return lcs;
}
