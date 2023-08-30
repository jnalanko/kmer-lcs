#include <sdsl/int_vector.hpp>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include <cmath>
#include <unordered_map>

char get_keyab_from_key(char key, const std::unordered_map<char,char>& alphabet_map) {
    auto it = alphabet_map.find(key);
    if (it != alphabet_map.end()) {
        return it->second;
    }
    return 0;  // Return a default value if key is not found
}

struct alphabet_str{
    vector<char> super_alphabet;
    unordered_map<char,char> super_alphabet_map;
};

alphabet_str get_super_alphabet(const vector<char>& alphabet, const char base){ // base = base of the original alphabet
    alphabet_str result;
    char excluded_v = (1ULL << (static_cast<char>(pow(3, base/2))-1));
    vector<char> super_alphabet;
    unordered_map<char, char> super_alphabet_map;
    super_alphabet_map.reserve(alphabet.size()*(alphabet.size()-1)+1);
    //remove chars w $
    for(char c2:alphabet){
        for (char c1:alphabet){
            char  super_char = ( c1 << (3 * base)) | c2;
            super_alphabet_map[super_char] = super_alphabet.size();
            super_alphabet.push_back(super_char);
            if (super_char < excluded_v){break;} // avoid char ending with $
        }
    }
    result.super_alphabet = super_alphabet;
    result.super_alphabet_map = super_alphabet_map;
    return result;
}

bool compare_2_char(char compare_res){
    compare_res = compare_res & 0b111;
    return compare_res == 0;
}

struct concat_str{
    vector<char> concat_v;
    sdsl::bit_vector concat_b;
};

concat_str get_super_concat(char base, const vector<char>& concat_vector, const sdsl::bit_vector& concat_bits, const vector<uint64_t>& C_array, const vector<char>& super_alphabet,const unordered_map<char,char>& super_alphabet_map, const vector<char>& alphabet){
    concat_str result;
    vector<char> res_concat_v;
    vector<bool> bv;

    // Scan concat_bv once and add a pointer for every element in the alphabet
    vector<uint64_t> select_v = {0};
    uint64_t i=1, subset_counter = 0; // In the first pos there is a 1 = {} = $$...
    for (uint64_t k = 0; k < C_array.size(); k++){
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

    vector<uint64_t> C_copy = {0};
    C_copy.insert(C_copy.end(), C_array.begin(), C_array.end());
    uint64_t subsets = 0;
    for (uint64_t x = 0; x < concat_vector.size(); x++){
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
        uint64_t x_subset = C_copy[curr_char];

        // Start of the subset of curr_char
        // INDEX IN BITVECTOR
        uint64_t x_subset_start_bv = select_v[curr_char] + 1; // later +1 as we want zeros
        // INDEX CONCAT VECTOR
        uint64_t x_subset_start_v = x_subset_start_bv - x_subset;//rank 0
        // Scan the subset and add each character found to curr_char
        // (2) Check you are not appending an empty set
        while (x_subset_start_bv < concat_bits.size() && !concat_bits[x_subset_start_bv]){
            char key = (alphabet[curr_char] << (3*(base/2)))| alphabet[concat_vector[x_subset_start_v]] ;
            char key_ab = get_keyab_from_key( key, super_alphabet_map);

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


sdsl::int_vector<> lcs_superalphabet_algorithm(const sbwt::plain_matrix_sbwt_t& SBWT, char chars_per_super_alpha_char){
    cerr << "Starting Superalphabet" << endl;
    char t = 2;
    const sdsl::bit_vector& A_bits = SBWT.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = SBWT.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = SBWT.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = SBWT.get_subset_rank_structure().T_bits;
    int64_t k = SBWT.get_k();
    uint64_t n_nodes = SBWT.number_of_subsets();
    vector<uint64_t> C_array(4);

    vector<char> last = {0};
    C_array[0] = last.size();
    for(uint64_t i = 0; i < n_nodes; i++) if(A_bits[i]) last.push_back(1);
    C_array[1] = last.size();
    for(uint64_t i = 0; i < n_nodes; i++) if(C_bits[i]) last.push_back(2);
    C_array[2] = last.size();
    for(uint64_t i = 0; i < n_nodes; i++) if(G_bits[i]) last.push_back(3);
    C_array[3] = last.size();
    for(uint64_t i = 0; i < n_nodes; i++) if(T_bits[i]) last.push_back(4);
    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }
    sdsl::bit_vector mismatch_found_marks(n_nodes, 0);
    mismatch_found_marks[0]=1;
    sdsl::int_vector<> lcs(n_nodes, 0, 64 - __builtin_clzll(k-1)); // Enough bits per element to store values from 0 to k-1

    // Propagate the labels one step forward in the graph
    vector<char> propagated(n_nodes, 0);
    uint64_t A_ptr = C_array[0];
    uint64_t C_ptr = C_array[1];
    uint64_t G_ptr = C_array[2];
    uint64_t T_ptr = C_array[3];
    for(uint64_t i = 0; i < n_nodes-1; i++){
        if(mismatch_found_marks[i+1] == 0 &&  last[i] != last[i+1]){ // last[i-1] has been modified already
            mismatch_found_marks[i+1] = 1;
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
    for(uint64_t i = 0; i < n_nodes; i++){
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


    vector<char> alphabet = {0,1,2,3,4};//$,A,C,G,T

    // Create a super alphabet
    auto[super_2_alphabet, super_2_alphabet_map] = get_super_alphabet(alphabet, t/2);

    // Create a 2-C_array and check the 2nd char
    vector<uint64_t> super_2_C_array(super_2_alphabet.size(), 0);
    for ( size_t i = 1; i < n_nodes; i++ ) { // start from 1 as later i-1 and $$ already ok
        // mismatch in the 1st char
        char key = (propagated[i] << 3) | last[i];
        if (mismatch_found_marks[i]){
            char c = super_2_alphabet_map[key];
            while(super_2_C_array[c]==0){
                super_2_C_array[c]=i;
                c--;
            }
            // mismatch in the 2nd char
        } else if (!mismatch_found_marks[i] && (propagated[i] ^ propagated[i-1])!=0){
            mismatch_found_marks[i] = 1;
            lcs[i] = 1;
            char c = super_2_alphabet_map[key];
            while(super_2_C_array[c]==0){
                super_2_C_array[c]=i;
                c--;
            }

        } // else no mismatch yet
        // Create the first 2 columns according to the new super_2_C_array

        last[i]= super_2_alphabet_map[key];//last[i] = super_2_alphabet.size()-1;//get_keyab_from_key(key, super_2_alphabet);
    }
    super_2_C_array.erase(super_2_C_array.begin());// remove the first 0


    // CONCAT 2-super alphabet
    auto[super_concat_v, super_concat_b] = get_super_concat(2,concat_v, concat_b, C_array, super_2_alphabet, super_2_alphabet_map, alphabet );

// -------------------------------------------------------------------
    // Dump k-mers with 2-super_alphabet
    // Propagate the labels one step forward in the graph

    for (uint64_t round = t; round < (k+t-1); round += t) { // works w chars_per_super_alpha_char = 2
        // Propagate the labels one step forward in the graph
        cerr << "Round " << to_string(round) << endl;

        vector <uint64_t> super_C_copy = {0};
        super_C_copy.insert(super_C_copy.end(), super_2_C_array.begin(), super_2_C_array.end());
        uint64_t subsets = 0;
        for (uint64_t i = 0; i < super_concat_v.size(); i++) {
            if (super_concat_b[i + subsets]) {
                subsets++;
                if (super_concat_b[i + subsets]) { // empty subset
                    i--;
                    continue;
                }
            }
            propagated[super_C_copy[super_concat_v[i]]++] = last[subsets - 1];
        }
        if (chars_per_super_alpha_char > t) {
            t *=2;
            break;}
        for (uint64_t j = 1; j < n_nodes; j++) {
            if (!mismatch_found_marks[j]) {
                char pi = super_2_alphabet[propagated[j]];
                char pi_ = super_2_alphabet[propagated[j - 1]];
                // Compare pi and pi-1 with xor and mask the result to see only the last 3 bits
                char compare_res = pi ^ pi_;
                if (compare_res != 0) {
                    mismatch_found_marks[j] = 1;
                    lcs[j] = round + compare_2_char(compare_res);
                }
            }
        }
        last = propagated;
    }

    if (chars_per_super_alpha_char == 2) {
        cerr << "Completed superalphabet-2"<< endl;
        return lcs;}
////--------------------------------------------------- t = 4
    //Create a 4-superalphabet
    auto[super_4_alphabet, super_4_alphabet_map] = get_super_alphabet(super_2_alphabet, t/2);

    // Create a 4-super_C_array
    vector<int64_t> super_4_C_array(super_4_alphabet.size(), 0);

    for ( size_t i = 1; i < n_nodes; i++ ) { // start from 1 as later i-1 and $$ already ok
        char pi = super_2_alphabet[propagated[i]];
        char li = super_2_alphabet[last[i]];
        char key = ( pi << (3*t/2)) | li;

        if (mismatch_found_marks[i]){
            char c = super_4_alphabet_map[key];
            while(super_4_C_array[c]==0){
                super_4_C_array[c]=i;
                c--;
            }
            // mismatch in the 3rd or 4th char
        } else if (!mismatch_found_marks[i]){
            // Compare pi and pi-1 with xor and mask the result to see only the last 3 bits
            char pi_ = super_2_alphabet[propagated[i-1]];
            char compare_res = pi ^ pi_;
            if (compare_res != 0) {
                mismatch_found_marks[i] = 1;
                lcs[i] = 2 + compare_2_char(compare_res);
                char c = super_4_alphabet_map[key];
                while(super_4_C_array[c]==0){
                    super_4_C_array[c]=i;
                    c--;
                }
            }
        }
        // Create the first 4 columns according to the new super_C_array
        last[i]= super_4_alphabet_map[key];
    }
    super_4_C_array.erase(super_4_C_array.begin());// remove the first 0

    int64_t c = super_4_C_array.size();
    while( super_4_C_array[c] == 0){
        super_4_C_array[c] = n_nodes;
        c--;
    }

    // CONCAT 4-super alphabet
    auto concat_4 = get_super_concat(4,super_concat_v, super_concat_b, super_2_C_array, super_4_alphabet,super_4_alphabet_map, super_2_alphabet );
    super_concat_v = concat_4.concat_v;
    super_concat_b = concat_4.concat_b;

////-------------
////     Dump k-mers with 4-super_alphabet
////     Propagate the labels one step forward in the graph
    for (uint64_t round = t; round < (k+t-1); round += t) { // works w chars_per_super_alpha_char = 2
        // Propagate the labels one step forward in the graph
        cerr << "Round_4 " << to_string(round) << endl;
        vector <uint64_t> super_C_copy = {0};
        super_C_copy.insert(super_C_copy.end(), super_4_C_array.begin(), super_4_C_array.end());
        uint64_t subsets = 0;
        for (uint64_t i = 0; i < super_concat_v.size(); i++) {
            if (super_concat_b[i + subsets]) {
                subsets++;
                if (super_concat_b[i + subsets]) { // empty subset
                    i--;
                    continue;
                }
            }
            propagated[super_C_copy[super_concat_v[i]]++] = last[subsets - 1];
        }
//        if (chars_per_super_alpha_char > t) {
//            t *=2;
//            break;}
        for (uint64_t j = 1; j < n_nodes; j++) {
            if (!mismatch_found_marks[j]) {
                char pi = super_4_alphabet[propagated[j]];
                char pi_ = super_4_alphabet[propagated[j-1]];
                // Compare pi and pi-1 with xor and mask the result to see only the last 3 bits
                char compare_res = pi ^ pi_;
                if (compare_res != 0) {
                    mismatch_found_marks[j] = 1;
                    // Check only c2 (last 6 bits)
                    char compare_res_12 = compare_res & 0b111111;
                    if (compare_res_12 != 0) { // the last 6 bits are not the same
                        lcs[j] = round + compare_2_char(compare_res_12);
                    } else { // No mismatch in the last 6 bits -> mismatch in the previous
                        lcs[j] = round + 2 + compare_2_char(compare_res >> 6);
                    }
                }
            }
        }
        last = propagated;
    }

    return lcs;
}

