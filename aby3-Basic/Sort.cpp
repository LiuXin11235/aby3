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

    // if(pIdx == 0){
    //     debug_info("bc sorts - num_interval = " + std::to_string(num_interval) + ", total_size = " + std::to_string(total_size) + ", valid_size = " + std::to_string(valid_size));
    // }

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

            // if(pIdx == 0) debug_info("current after size = " + std::to_string((after - total_current)) + " | after = " + std::to_string(after) + " | total_current = " + std::to_string(total_current));
            // if(pIdx == 0) debug_info("i = " + std::to_string(i) + " | low = " + std::to_string(low) + " | high = " + std::to_string(high) + " | current = " + std::to_string(current));

            if((after - total_current) > valid_size) break;

            // tiling the indices.
            std::vector<size_t> ranging_indices(high - low, 0);
            std::iota(ranging_indices.begin(), ranging_indices.end(), low);

            std::string ranging_indices_str = "";
            for(size_t j=0; j<(high - low); j++){
                ranging_indices_str += std::to_string(ranging_indices[j]) + " ";
            }    

            vector_tile<size_t>(ranging_indices, (high - low), index_tile, current - total_current);
            vector_repeat<size_t>(ranging_indices, (high - low), index_repeat, current - total_current);
            current = after;
        }

        if(pIdx == 0){
            std::string ranging_indices_str = "";
            for(size_t i=0; i<index_tile.size(); i++){
                ranging_indices_str += std::to_string(index_tile[i]) + " ";
            }
        }

        // if(pIdx == 0) debug_info("bc sorts - before tiling & repeating data construction.");

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

        bool_cipher_lt(pIdx, data_tile, data_repeat, data_comp, enc, eval, runtime);
        i64Matrix plain_comp(current - total_current, 1);
        enc.revealAll(runtime, data_comp, plain_comp).get();

        // compute the location after sorting and construct the data.
        current = total_current;
        size_t next_start = start;

        for(size_t i=0; i<(num_interval - start); i++){
            size_t low = lows[start + i], high = highs[start + i];
            size_t after = current + (high - low) * (high - low);

            if((after - total_current) > valid_size){
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


int quick_sort_different(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size){
    size_t BITSIZE = data[0].bitCount();
    size_t len = data.size();

    std::vector<size_t> lows = {0};
    std::vector<size_t> highs = {len};
    size_t num_interval = 1;

    std::vector<size_t> bc_lows;
    std::vector<size_t> bc_highs;

    if(min_size < 4) min_size = 4;

    while(num_interval > 0){
        std::vector<size_t> mids(lows.size());
        for(size_t i=0; i<lows.size(); i++){
            mids[i] = ((highs[i] - 1 - lows[i]) >> 1) + lows[i];
        }
        
        // Select out the pivot candidates for each 
        std::vector<sbMatrix> pivot_candidates(num_interval * 3);
        for(size_t i=0; i<num_interval; i++){
            pivot_candidates[i * 3] = data[lows[i]];
            pivot_candidates[i * 3 + 1] = data[mids[i]];
            pivot_candidates[i * 3 + 2] = data[highs[i] - 1];
        }

        sbMatrix tile_candidate_pivots(num_interval * 3 * 3, BITSIZE);
        sbMatrix repeat_candidate_pivots(num_interval * 3 * 3, BITSIZE);

        for(size_t i=0; i < num_interval; i++){
            size_t start_ind = i * 3 * 3;
            for(size_t j=0; j<3; j++){
                for(size_t k=0; k<3; k++){
                    tile_candidate_pivots.mShares[0](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + k].mShares[0](0, 0);
                    tile_candidate_pivots.mShares[1](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + k].mShares[1](0, 0);
                    repeat_candidate_pivots.mShares[0](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + j].mShares[0](0, 0);
                    repeat_candidate_pivots.mShares[1](start_ind + j*3 + k, 0) = pivot_candidates[i * 3 + j].mShares[1](0, 0);
                }
            }
        }

        // compare the pivots in each interval through one vectorized lt and then reveal it.
        sbMatrix comp_candidate_pivots(num_interval * 3 * 3, 1);
        bool_cipher_lt(pIdx, tile_candidate_pivots, repeat_candidate_pivots, comp_candidate_pivots, enc, eval, runtime);
        i64Matrix plain_comp_candidates(num_interval * 3 * 3, 1);
        enc.revealAll(runtime, comp_candidate_pivots, plain_comp_candidates).get();

        // compute the indices of each interval pivots.
        std::vector<std::vector<int>> comp_sum(num_interval);
        for(size_t i=0; i<num_interval; i++){
            comp_sum[i].resize(3);
            std::fill(comp_sum[i].begin(), comp_sum[i].begin() + 3, 0);
            for(size_t j=0; j<3; j++){
                for(size_t k=0; k<3; k++) comp_sum[i][j] += (int) plain_comp_candidates(i * 3*3 + j * 3 + k, 0);
            }
        }

        // compute the indices of each interval pivots.
        std::vector<size_t> which_pivots = argwhere(comp_sum, 1);
        
        // swicth the elements of each interval, put the pivot in the first.
        std::vector<size_t> switch_mid = argwhere(which_pivots, 1);
        for(size_t i=0; i<switch_mid.size(); i++){
            size_t si = switch_mid[i];
            data[lows[si]] = pivot_candidates[si * 3 + 1];
            data[mids[si]] = pivot_candidates[si * 3];
        }

        std::vector<size_t> switch_high = argwhere(which_pivots, 2);
        for(size_t i=0; i<switch_high.size(); i++){
            size_t si = switch_high[i];
            data[lows[si]] = pivot_candidates[si * 3 + 2];
            data[highs[si] - 1] = pivot_candidates[si * 3];
        }

        // vectorized compare the elements with the midpoint of each interval.
        size_t total_cmp_size = std::accumulate(highs.begin(), highs.end(), 0) - std::accumulate(lows.begin(), lows.end(), 0) - num_interval;
        std::vector<size_t> to_cmp_with_low_pos(total_cmp_size, 0);
        std::vector<size_t> to_cmp_low_pos(total_cmp_size, 0);

        size_t current = 0;
        for(size_t i=0; i<num_interval; i++){
            size_t low = lows[i], high = highs[i];
            size_t after = current + high - low - 1;

            for(size_t j=current; j<after; j++){
                to_cmp_with_low_pos[j] = low + 1 + j - current;
                to_cmp_low_pos[j] = low;
            }
            current = after;
        }


        sbMatrix to_cmp_low(total_cmp_size, BITSIZE);
        sbMatrix to_cmp_with_low(total_cmp_size, BITSIZE);
        sbMatrix tmp_cmp(total_cmp_size, 1);
        i64Matrix plain_cmp(total_cmp_size, 1);

        for(size_t i=0; i<total_cmp_size; i++){
            to_cmp_low.mShares[0](i, 0) = data[to_cmp_low_pos[i]].mShares[0](0, 0);
            to_cmp_low.mShares[1](i, 0) = data[to_cmp_low_pos[i]].mShares[1](0, 0);
            to_cmp_with_low.mShares[0](i, 0) = data[to_cmp_with_low_pos[i]].mShares[0](0, 0);
            to_cmp_with_low.mShares[1](i, 0) = data[to_cmp_with_low_pos[i]].mShares[1](0, 0);
        }
        bool_cipher_lt(pIdx, to_cmp_low, to_cmp_with_low, tmp_cmp, enc, eval, runtime);
        enc.revealAll(runtime, tmp_cmp, plain_cmp).get();

        // update the intervals.
        std::vector<size_t> new_lows;
        std::vector<size_t> new_highs;
        std::vector<size_t> data_pos(current + num_interval, 0);
        std::vector<size_t> origion_pos(current + num_interval, 0);

        // construct the corresponsing positions.
        current = 0;
        size_t new_current = 0;
        for(size_t i=0; i<num_interval; i++){
            size_t low = lows[i], high = highs[i];
            size_t after = current + high - low - 1;
            size_t new_after = new_current + high - low;

            size_t large_num = 0;
            for(size_t j=current; j<after; j++){
                large_num += (size_t) plain_cmp(j, 0);
            }
            size_t low_pos = high - 1 - large_num;
            for(size_t j=new_current; j<new_after; j++) data_pos[j] = low + j - new_current;

            aby3::i64Matrix cmp_interval = plain_cmp.block(current, 0, after - current, 1);

            std::vector<size_t> tmp_lower_args = argwhere(cmp_interval, 0);
            std::vector<size_t> tmp_higher_args = argwhere(cmp_interval, 1);

            for(size_t j=new_current; j<new_after; j++){
                if(j < new_after - large_num - 1){
                    origion_pos[j] = low + 1 + tmp_lower_args[j - new_current];
                }
                else if(j > new_after - large_num - 1){
                    origion_pos[j] = low + 1 + tmp_higher_args[j - (new_after - large_num)];
                }
                else{
                    origion_pos[j] = low;
                }
            }

            // construct the following interval end-points.
            if((low_pos - low) > min_size){
                new_lows.push_back(low);
                new_highs.push_back(low_pos);
            }
            else{
                if((low_pos - low) > 1){
                    bc_lows.push_back(low);
                    bc_highs.push_back(low_pos);
                }
            }

            if((high - low_pos - 1) > min_size){
                new_lows.push_back(low_pos + 1);
                new_highs.push_back(high);
            }
            else{
                if((high - low_pos - 1) > 1) {
                    bc_lows.push_back(low_pos + 1);
                    bc_highs.push_back(high);
                }
            }

            current = after;
            new_current = new_after;
        }
        
        std::vector<sbMatrix> tmp_data = data;
        for(size_t i=0; i<current + num_interval; i++){
            data[data_pos[i]] = tmp_data[origion_pos[i]];
        }

        lows = new_lows;
        highs = new_highs;
        num_interval = lows.size();
    }

    // if(pIdx == 0) debug_info("before bc_sort!");
    // if(pIdx == 0){
    //     for(size_t i=0; i<bc_lows.size(); i++){
    //         debug_info("bc_lows[" + std::to_string(i) + "] = " + std::to_string(bc_lows[i]) + ", bc_highs[" + std::to_string(i) + "] = " + std::to_string(bc_highs[i]) + " interval_size = " + std::to_string(bc_highs[i] - bc_lows[i]));
    //     }
    // }

    bc_sort_different(data, bc_lows, bc_highs, pIdx, enc, eval, runtime, 1048576);
    return 0;
}

int quick_sort(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size){

    tag_append(pIdx, data);
    quick_sort_different(data, pIdx, enc, eval, runtime, min_size);
    size_t tag_size = std::ceil(std::log2(data.size()));
    tag_remove(pIdx, tag_size, data);

    return 0;
}