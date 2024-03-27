#include "DataGenerator.h"

std::vector<aby3::si64> generate_vector_si64(size_t len, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime){
    std::vector<aby3::si64> target_vec(len);
    for(size_t k=0; k<len; k+=ENCRYPTION_LIMIT){
        size_t start = k, end = std::min(k+ENCRYPTION_LIMIT, len);
        size_t block_size = end - start;
        aby3::si64Matrix target_mat(block_size, 1);
        init_zeros(pIdx, enc, runtime, target_mat, block_size);
        for(size_t i=start; i<end-start; i++) target_vec[i] = target_mat(i-start, 0);
    }
    return target_vec;
}
std::vector<int> generate_vector_int(size_t len){
    std::vector<int> target_vec(len);
    for(size_t i=0; i<len; i++) target_vec[i] = i;
    return target_vec;
}