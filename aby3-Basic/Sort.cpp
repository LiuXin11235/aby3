#include "Sort.h"
#include <numeric>

using namespace oc;
using namespace aby3;

int bc_sort_different(std::vector<aby3::sbMatrix> &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size){

    size_t BITSIZE = data[0].bitCount();

    size_t num_interval = lows.size();
    size_t total_size = 0;
    for(size_t i=0; i<num_interval; i++){
        total_size += (highs[i] - lows[i]) * (highs[i] - lows[i]);
    }
    size_t valid_size = std::min(total_size, max_size);

    // Prepare the tiling & repeating indices.
    std::vector<size_t> index_tile(valid_size, 0);
    std::vector<size_t> index_repeat(valid_size, 0);

    // Setup the recorders.
    size_t start = 0;
    size_t total_current = 0;
    size_t pos_current = 0;

    // Prepare the data location exchanging vectors.
    size_t pos_num = std::accumulate(highs.begin(), highs.end(), 0) - std::accumulate(lows.begin(), lows.end(), 0);

    std::vector<size_t> data_pos(pos_num, 0);
    std::vector<size_t> origion_pos(pos_num, 0);

    while(start < num_interval){
        size_t current = total_current;
        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            if(after > valid_size) break;

            // tiling the indices.
            std::vector<size_t> ranging_indices(high - low, low);
            std::iota(ranging_indices.begin(), ranging_indices.end(), low);

            vector_tile<size_t>(ranging_indices, (high - low), index_tile, current - total_current);
            vector_repeat<size_t>(ranging_indices, (high - low), index_repeat, current - total_current);
            current = after;
        }

        // vectorized comparison.
        sbMatrix data_tile(current - total_current, BITSIZE);
        sbMatrix data_repeat(current - total_current, BITSIZE);
        sbMatrix data_comp(current - total_current, 1);
        for(size_t i=0; i<(current - total_current); i++){
            data_tile.mShares[0](i, 0) = data[index_tile[i]].mShares[0](0, 0);
            data_tile.mShares[1](i, 0) = data[index_tile[i]].mShares[1](0, 0);
            data_repeat.mShares[0](i, 0) = data[index_repeat[i]].mShares[0](0, 0);
            data_repeat.mShares[1](i, 0) = data[index_repeat[i]].mShares[1](0, 0);
        }

        i64Matrix debug_mat(current - total_current, 1);
        enc.revealAll(runtime, data_tile, debug_mat).get();
        bool_cipher_lt(pIdx, data_tile, data_repeat, data_comp, enc, eval, runtime);
        i64Matrix plain_comp(current - total_current, 1);
        enc.revealAll(runtime, data_comp, plain_comp).get();

        // compute the location after sorting and construct the data.
        current = total_current;
        size_t next_start = start;

        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            if(after > valid_size){
                next_start = start + i;
                break;
            }

            next_start += 1;
            size_t pos_after = pos_current + high - low;
            std::vector<size_t> rela_pos(high - low, low);
            for(size_t i=0; i<high - low; i++){
                for(size_t j=0; j<high - low; j++){
                    rela_pos[i] += (size_t) plain_comp(current - total_current + i * (high - low) + j, 0);
                }
            }

            for(size_t i=pos_current; i<pos_after; i++){
                data_pos[i] = rela_pos[i - pos_current];
                origion_pos[i] = i - pos_current + low;
            }

            pos_current = pos_after;
            current = after;
        }

        start = next_start;
        total_current = current;
        
    }

    // sort the data.
    std::vector<sbMatrix> tmp_data = data;
    for(size_t i=0; i<pos_num; i++){
        data[data_pos[i]] = tmp_data[origion_pos[i]];
    }   

    return 0;
}