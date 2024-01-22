#include "Basics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <random>

#include <aby3/sh3/Sh3BinaryEvaluator.h>
#include <aby3/Circuit/CircuitLibrary.h>

#include "../aby3-RTR/debug.h"

// #define DEBUG_BASIC
static std::string PARTY_FILE = "/root/aby3/party-";

using namespace oc;
using namespace aby3;

void bool_cipher_lt(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_lt(bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedB);
    binEng.setInput(1, sharedA);    
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, 1);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_eq(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_eq(bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, 1);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_or(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_bitwiseOr(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}


void bool_cipher_add(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_add(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_and(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_bitwiseAnd(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_init_false(int pIdx, aby3::sbMatrix &res){
    int bitSize = res.bitCount();
    int i64Size = res.i64Size();

    for(int i=0; i<i64Size; i++){
        switch(pIdx){
            case 0:
                res.mShares[0](i, 0) = 1;
                res.mShares[1](i, 0) = 0;
                break;
            case 1:
                res.mShares[0](i, 0) = 1;
                res.mShares[1](i, 0) = 1;
                break;
            case 2:
                res.mShares[0](i, 0) = 0;
                res.mShares[1](i, 0) = 1;
                break;
            default:
                throw std::runtime_error("bool_init_false: pIdx out of range.");
        }
    }
    return;
}

void bool_init_true(int pIdx, aby3::sbMatrix &res){
    int bitSize = res.bitCount();
    int i64Size = res.i64Size();

    for(int i=0; i<i64Size; i++){
        switch(pIdx){
            case 0:
                res.mShares[0](i, 0) = 0;
                res.mShares[1](i, 0) = 0;
                break;
            case 1:
                res.mShares[0](i, 0) = 1;
                res.mShares[1](i, 0) = 0;
                break;
            case 2:
                res.mShares[0](i, 0) = 0;
                res.mShares[1](i, 0) = 1;
                break;
            default:
                throw std::runtime_error("bool_init_true: pIdx out of range.");
        }
    }
    return;
}

void bool_init_false(int pIdx, boolShare &res){
    switch(pIdx){
        case 0:
            res.bshares[0] = 1;
            res.bshares[1] = 0;
            break;
        case 1:
            res.bshares[0] = 1;
            res.bshares[1] = 1;
            break;
        case 2:
            res.bshares[0] = 0;
            res.bshares[1] = 1;
            break;
        default:
            throw std::runtime_error("bool_init_false: pIdx out of range.");
    }
    return;
}

void bool_init_true(int pIdx, boolShare &res){
    switch(pIdx){
        case 0:
            res.bshares[0] = 0;
            res.bshares[1] = 0;
            break;
        case 1:
            res.bshares[0] = 1;
            res.bshares[1] = 0;
            break;
        case 2:
            res.bshares[0] = 0;
            res.bshares[1] = 1;
            break;
        default:
            throw std::runtime_error("bool_init_true: pIdx out of range.");
    }
    return;
}


void bool_shift_and_left(int pIdx, aby3::sbMatrix &sharedA, size_t shift_len, aby3::sbMatrix &res_shift, aby3::sbMatrix &res_left){
    // shift the input.
    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();
    std::cout << "bitSize: " << bitSize << std::endl;
    int left_size = bitSize - shift_len;

    // right shift each input to shift_len bits.
    for(size_t i=0; i<i64Size; i++){
        res_shift.mShares[0](i, 0) = sharedA.mShares[0](i, 0) >> shift_len;
        res_shift.mShares[1](i, 0) = sharedA.mShares[1](i, 0) >> shift_len;
    }

    // left shift each input to bitSize - shift_len bits.
    for(size_t i=0; i<i64Size; i++){
        res_left.mShares[0](i, 0) = sharedA.mShares[0](i, 0) & ((1 << shift_len) - 1);
        res_left.mShares[1](i, 0) = sharedA.mShares[1](i, 0) & ((1 << shift_len) - 1);
    }

    // fulfill the leading bits to logical 0.
    switch(pIdx) {
        case 0:
            for(size_t i=0; i<i64Size; i++){
                res_shift.mShares[0](i, 0) = res_shift.mShares[0](i, 0) | (~0 << shift_len);
                res_shift.mShares[1](i, 0) = res_shift.mShares[1](i, 0) & ~(~0 << shift_len);
                res_left.mShares[0](i, 0) = res_left.mShares[0](i, 0) | (~0 << left_size);
                res_left.mShares[1](i, 0) = res_left.mShares[1](i, 0) & ~(~0 << left_size);
            }
            break;
        case 1:
            for(size_t i=0; i<i64Size; i++){
                res_shift.mShares[0](i, 0) = res_shift.mShares[0](i, 0) | (~0 << shift_len);
                res_shift.mShares[1](i, 0) = res_shift.mShares[1](i, 0) | (~0 << shift_len);
                res_left.mShares[0](i, 0) = res_left.mShares[0](i, 0) | (~0 << left_size);
                res_left.mShares[1](i, 0) = res_left.mShares[1](i, 0) | (~0 << left_size);
            }
            break;
        case 2:
            for(size_t i=0; i<i64Size; i++){
                res_shift.mShares[0](i, 0) = res_shift.mShares[0](i, 0) & ~(~0 << shift_len);
                res_shift.mShares[1](i, 0) = res_shift.mShares[1](i, 0) | (~0 << shift_len);
                res_left.mShares[0](i, 0) = res_left.mShares[0](i, 0) & ~(~0 << left_size);
                res_left.mShares[1](i, 0) = res_left.mShares[1](i, 0) | (~0 << left_size);
            }
            break;
        default:
            throw std::runtime_error("bool_shift_and_left: pIdx out of range.");
    }
}

aby3::i64Matrix back2plain(int pIdx, aby3::sbMatrix &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    // send the shares to the prevParty.
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.mShares[0]);

    // rece from the nextParty.
    aby3::i64Matrix other_share(cipher_val.i64Size(), 1);
    runtime.mComm.mNext.recv(other_share.data(), other_share.size());

    // get the shares.
    for(size_t i=0; i<cipher_val.i64Size(); i++){
        other_share(i, 0) = other_share(i, 0) ^ cipher_val.mShares[1](i, 0) ^ cipher_val.mShares[0](i, 0);
    }
    return other_share;
}


bool back2plain(int pIdx, boolShare &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.bshares[0]);
    bool other_share;
    runtime.mComm.mNext.recv(other_share);
    other_share = other_share ^ cipher_val.bshares[1] ^ cipher_val.bshares[0];
    return other_share;
}


aby3::i64 back2plain(int pIdx, boolIndex &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
#ifdef DEBUG_BASIC
    std::string debug_party_file = PARTY_FILE + std::to_string(pIdx);
    std::ofstream ofs(debug_party_file, std::ios::app);
    debug_info("share0: " + std::to_string(cipher_val.indexShares[0]), ofs);
    debug_info("share1: " + std::to_string(cipher_val.indexShares[1]), ofs);
#endif
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.indexShares[0]);
    aby3::i64 other_share = 0;
    runtime.mComm.mNext.recv(other_share);

#ifdef DEBUG_BASIC
    debug_info("other share: " + std::to_string(other_share), ofs);
#endif
    other_share = other_share ^ cipher_val.indexShares[1] ^ cipher_val.indexShares[0];
    return other_share;
}


void get_permutation(size_t len, std::vector<size_t> &permutation, block &seed){
    permutation.resize(len);
    for(size_t i = 0; i < len; i++){
            permutation[i] = i;
    }
    PRNG prng(seed);    
    std::random_shuffle(permutation.begin(), permutation.end(), prng);
    return;
}

void get_inverse_permutation(std::vector<size_t> &permutation, std::vector<size_t> &inverse_permutation){
    size_t len = permutation.size();
    inverse_permutation.resize(len);
    for(size_t i = 0; i < len; i++){
        inverse_permutation[permutation[i]] = i;
    }
    return;
}

void combine_permutation(std::vector<std::vector<size_t>> &permutation_list, std::vector<size_t> &final_permutation){
    size_t len = permutation_list[0].size();
    size_t permutes = permutation_list.size();
    final_permutation.resize(len);
    for(size_t i = 0; i < len; i++){
        final_permutation[i] = i;
    }
    for(size_t i = 0; i < permutes; i++){
        std::vector<size_t> tmp_inverse = permutation_list[permutes - i - 1];
        get_inverse_permutation(permutation_list[permutes - i - 1], tmp_inverse);
        plain_permutate(tmp_inverse, final_permutation);
    }
    return;
}

void get_random_mask(int pIdx, i64Matrix &res, block &seed){
    size_t len = res.rows();
    PRNG prng(seed);
    for(size_t i=0; i<len; i++) res(i, 0) = (i64) prng.get<int64_t>();
    return;
}