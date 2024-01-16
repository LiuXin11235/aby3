#include "Shuffle.h"

using namespace oc;
using namespace aby3;

// #define DEBUG_SHUFFLE
static std::string PARTY_FILE = "/root/aby3/party-";


/**
 * The shuffle protocol accoring to https://dl.acm.org/doi/10.1145/3460120.3484560 (advanced version).
*/
int efficient_shuffle(std::vector<sbMatrix> &T, int pIdx, std::vector<sbMatrix> &Tres, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.size();
    size_t unit_len = T[0].i64Size();

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    std::string debug_party_file = PARTY_FILE + std::to_string(pIdx) + ".txt";
    std::ofstream ofs(debug_party_file, std::ios_base::app);

#ifdef DEBUG_SHUFFLE
    std::vector<size_t> other_permutation(len);
    runtime.mComm.mPrev.asyncSendCopy(next_permutation.data(), next_permutation.size());
    runtime.mComm.mNext.recv(other_permutation.data(), other_permutation.size());

    debug_output_vector(prev_permutation, ofs);
    debug_output_vector(next_permutation, ofs);
    debug_output_vector(other_permutation, ofs);

    std::vector<size_t> final_permutation;  
    std::vector<std::vector<size_t>> permutation_list;
    if(pIdx == 0){
        permutation_list = {next_permutation, prev_permutation, other_permutation};
    }
    if(pIdx == 1){
        permutation_list = {prev_permutation, other_permutation, next_permutation};
    }
    if(pIdx == 2){
        permutation_list = {other_permutation, next_permutation, prev_permutation};
    }
    // std::vector<std::vector<size_t>> permutation_list = {prev_permutation, next_permutation, other_permutation};
    combine_permutation(permutation_list, final_permutation);
    debug_output_vector(final_permutation, ofs);
#endif

    // 2 - generate the random masks Z.
    std::vector<i64Matrix> prev_maskZ(len);
    std::vector<i64Matrix> next_maskZ(len);
    for(size_t i = 0; i<len; i++){
        prev_maskZ[i].resize(unit_len, T[0].bitCount());
        next_maskZ[i].resize(unit_len, T[0].bitCount());
        get_random_mask(pIdx, prev_maskZ[i], prevSeed);
        get_random_mask(pIdx, next_maskZ[i], nextSeed);
    }

    // get the random permutation.
    if(pIdx == 0){
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        std::vector<i64Matrix> maskA(len);
        for(size_t i = 0; i<len; i++){
            maskB[i].resize(unit_len, T[0].bitCount());
            maskA[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskB[i], nextSeed);
            get_random_mask(pIdx, maskA[i], prevSeed);
        }

        // get the sharedX1.
        std::vector<i64Matrix> sharedX1(len);
        for(size_t i=0; i<len; i++){
            sharedX1[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                sharedX1[i](j, 0) = T[i].mShares[0](j) ^ T[i].mShares[1](j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX1);

        // get the sharedX2.
        std::vector<i64Matrix> sharedX2(len);
        for(size_t i=0; i<len; i++){
            sharedX2[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                sharedX2[i](j, 0) = sharedX1[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedX2);

        // send the sharedX2 to P1.
        i64Matrix sharedX2_matrix(len * unit_len, T[0].bitCount());
        for(size_t i=0; i<len; i++){
            for(size_t j=0; j<unit_len; j++){
                sharedX2_matrix(i * unit_len + j, 0) = sharedX2[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedX2_matrix.data(), sharedX2_matrix.size());

        // compute the final shares.
        for(size_t i=0; i<len; i++){
            Tres[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                Tres[i].mShares[1](j, 0) = maskA[i](j);
                Tres[i].mShares[0](j, 0) = maskB[i](j);
            }
        }
    }
    if(pIdx == 1){
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        for(size_t i = 0; i<len; i++){
            maskB[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskB[i], prevSeed);
        }

        // compute the sharedY1 and send to the next party.
        std::vector<i64Matrix> shared_Y1(len);
        for(size_t i=0; i<len; i++){
            shared_Y1[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, shared_Y1);

        i64Matrix sharedY1_matrix(len * unit_len, T[0].bitCount());
        for(size_t i=0; i<len; i++){
            for (size_t j = 0; j < unit_len; j++)
            {
                sharedY1_matrix(i*unit_len + j, 0) = shared_Y1[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedY1_matrix.data(), sharedY1_matrix.size());

        i64Matrix sharedX2_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(sharedX2_matrix.data(), sharedX2_matrix.size());

        // compute the sharedX3.
        std::vector<i64Matrix> sharedX3(len);
        for(size_t i=0; i<len; i++){
            sharedX3[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                sharedX3[i](j, 0) = sharedX2_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX3);

        // compute the masked C1.
        // std::vector<i64Matrix> maskedC1(len);
        i64Matrix maskedC1_matrix(len * unit_len, T[0].bitCount());
        for(size_t i=0; i<len; i++){
            for(size_t j=0; j<unit_len; j++){
                maskedC1_matrix(i*unit_len + j, 0) = sharedX3[i](j) ^ maskB[i](j);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(maskedC1_matrix.data(), maskedC1_matrix.size());

        i64Matrix maskedC2_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mNext.recv(maskedC2_matrix.data(), maskedC2_matrix.size());

        // compute the final shares.
        for(size_t i=0; i<len; i++){
            Tres[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                Tres[i].mShares[1](j, 0) = maskB[i](j);
                Tres[i].mShares[0](j, 0) = maskedC1_matrix(i*unit_len + j) ^ maskedC2_matrix(i*unit_len + j);
            }
        }

        
    }
    if(pIdx == 2){
        // pre-generate the randomness.
        std::vector<i64Matrix> maskA(len);
        for(size_t i=0; i<len; i++){
            maskA[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskA[i], nextSeed);
        }

        // receive the sharedY1 from the previous party.
        i64Matrix sharedY1_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(sharedY1_matrix.data(), sharedY1_matrix.size());

        // compute the sharedY2.
        std::vector<i64Matrix> sharedY2(len);
        for(size_t i=0; i<len; i++){
            sharedY2[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                sharedY2[i](j, 0) = sharedY1_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedY2);

        std::vector<i64Matrix> sharedY3(len);
        for(size_t i=0; i<len; i++){
            sharedY3[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                sharedY3[i](j, 0) = sharedY2[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedY3);


        // compute the masked C2.
        i64Matrix maskedC2_matrix(len * unit_len, T[0].bitCount());
        for(size_t i=0; i<len; i++){
            for(size_t j=0; j<unit_len; j++){
                maskedC2_matrix(i*unit_len + j, 0) = sharedY3[i](j) ^ maskA[i](j);
            }
        }
        runtime.mComm.mPrev.asyncSendCopy(maskedC2_matrix.data(), maskedC2_matrix.size());
        i64Matrix maskedC1_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(maskedC1_matrix.data(), maskedC1_matrix.size());

        std::vector<i64Matrix> maskC(len);
        for(size_t i=0; i<len; i++){
            maskC[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                maskC[i](j, 0) = maskedC1_matrix(i*unit_len + j) ^ maskedC2_matrix(i*unit_len + j);
            }
        }

        // compute the final shares.
        for(size_t i=0; i<len; i++){
            Tres[i].resize(unit_len, T[0].bitCount());
            for(size_t j=0; j<unit_len; j++){
                Tres[i].mShares[1](j, 0) = maskC[i](j);
                Tres[i].mShares[0](j, 0) = maskA[i](j);
            }
        }
    }
    return 0;
}