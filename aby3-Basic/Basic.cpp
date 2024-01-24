#include "Basics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <random>

#include <aby3/sh3/Sh3BinaryEvaluator.h>
#include <aby3/Circuit/CircuitLibrary.h>

#include "../aby3-RTR/debug.h"

#define DEBUG_BASIC
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
    res.resize(i64Size, bitSize);
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
    res.resize(i64Size, bitSize);
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


void bool_cipher_or(int pIdx, boolShare &sharedA, boolShare &sharedB, boolShare &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    bool cross1 = sharedA.bshares[0] && sharedB.bshares[0];
    bool cross2 = sharedA.bshares[0] && sharedB.bshares[1];
    bool cross3 = sharedA.bshares[1] && sharedB.bshares[0];
    bool share = (cross1 ^ cross2 ^ cross3);

    share = share ^ sharedA.bshares[0] ^ sharedB.bshares[0];

    runtime.mComm.mNext.asyncSendCopy(share);
    bool other_share;
    runtime.mComm.mPrev.recv(other_share);

    res.bshares = {share, other_share};

    return;
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

    return;
}

void bool_cipher_and(int pIdx, boolShare &sharedA, boolShare &sharedB, boolShare &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){

    bool cross1 = sharedA.bshares[0] && sharedB.bshares[0];
    bool cross2 = sharedA.bshares[0] && sharedB.bshares[1];
    bool cross3 = sharedA.bshares[1] && sharedB.bshares[0];
    bool share = (cross1 ^ cross2 ^ cross3);

    runtime.mComm.mNext.asyncSendCopy(share);
    bool other_share;
    runtime.mComm.mPrev.recv(other_share);

    res.bshares = {share, other_share};

    return;
}




void bool_cipher_not(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &res){

    int i64Size = sharedA.i64Size();

    switch(pIdx){
        case 0:
            for(size_t i=0; i<i64Size; i++){
                res.mShares[0](i, 0) = sharedA.mShares[0](i, 0);
                res.mShares[1](i, 0) = sharedA.mShares[1](i, 0);
            }
            break;
        case 1:
            for(size_t i=0; i<i64Size; i++){
                res.mShares[0](i, 0) = ~sharedA.mShares[0](i, 0);
                res.mShares[1](i, 0) = sharedA.mShares[1](i, 0);
            }
            break;
        case 2:
            for(size_t i=0; i<i64Size; i++){
                res.mShares[0](i, 0) = sharedA.mShares[0](i, 0);
                res.mShares[1](i, 0) = ~sharedA.mShares[1](i, 0);
            }
            break;
        default:
            throw std::runtime_error("bool_cipher_not: pIdx out of range.");
    }
    
    return;
}


void bool_cipher_not(int pIdx, std::vector<boolShare> &sharedA, std::vector<boolShare> &res){

    if(res.size() != sharedA.size()) res.resize(sharedA.size());
    size_t len = sharedA.size();
    switch(pIdx){
        case 0:
            for(size_t i=0; i<len; i++){
                res[i].bshares[0] = sharedA[i].bshares[0];
                res[i].bshares[1] = sharedA[i].bshares[1];
            }
            break;
        case 1:
            for(size_t i=0; i<len; i++){
                res[i].bshares[0] = !sharedA[i].bshares[0];
                res[i].bshares[1] = sharedA[i].bshares[1];
            }
            break;
        case 2:
            for(size_t i=0; i<len; i++){
                res[i].bshares[0] = sharedA[i].bshares[0];
                res[i].bshares[1] = !sharedA[i].bshares[1];
            }
            break;
        default:
            throw std::runtime_error("bool_cipher_not: pIdx out of range.");
    }
    return;
}

void bool_cipher_not(int pIdx, boolShare &sharedA, boolShare &res){
    switch(pIdx){
        case 0:
            res.bshares[0] = sharedA.bshares[0];
            res.bshares[1] = sharedA.bshares[1];
            break;
        case 1:
            res.bshares[0] = !sharedA.bshares[0];
            res.bshares[1] = sharedA.bshares[1];
            break;
        case 2:
            res.bshares[0] = sharedA.bshares[0];
            res.bshares[1] = !sharedA.bshares[1];
            break;
        default:
            throw std::runtime_error("bool_cipher_not: pIdx out of range.");
    }
    return;
}

void bool_cipher_dot(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    aby3::sbMatrix mul_res(i64Size, bitSize);

    auto cir = lib.int_int_bitwiseAnd(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        mul_res.resize(i64Size, bitSize);
        binEng.getOutput(0, mul_res);
    });
    dep.get();

    // std::cout << " after circuit eval" << std::endl;

    res.resize(1, bitSize);
    res.mShares[0](0, 0) = mul_res.mShares[0](0, 0);
    res.mShares[1](0, 0) = mul_res.mShares[1](0, 0);
    for (size_t i = 1; i < i64Size; i++)
    {
        res.mShares[0](0, 0) ^= mul_res.mShares[0](i, 0);
        res.mShares[1](0, 0) ^= mul_res.mShares[1](i, 0);
    }
}


void bool_get_first_zero_mask(int pIdx, std::vector<boolShare>& inputA, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    
    size_t len = inputA.size();
    std::vector<boolShare> tmp_mask(len);

    bool_cipher_not(pIdx, inputA, tmp_mask);

    boolShare flag = boolShare(0, 0);
    boolShare nFlag;

    // TODO: avoid the sequential computation?
    for(size_t i=0; i<len; i++){
        // flag = tmp_mask[i] & !flag;
        boolShare tmp;
        bool_cipher_not(pIdx, flag, nFlag);
        bool_cipher_and(pIdx, tmp_mask[i], nFlag, tmp, enc, eval, runtime);
        bool_cipher_or(pIdx, tmp, flag, flag, enc, eval, runtime);

        res.mShares[0](i, 0) = (int) tmp.bshares[0];
        res.mShares[1](i, 0) = (int) tmp.bshares[1];
    }
    return;
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
    runtime.mComm.mPrev.asyncSendCopy(cipher_val.indexShares[0]);
    aby3::i64 other_share = 0;
    runtime.mComm.mNext.recv(other_share);

    other_share = other_share ^ cipher_val.indexShares[1] ^ cipher_val.indexShares[0];
    return other_share;
}


std::vector<bool> back2plain(int pIdx, std::vector<boolShare> &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    std::vector<char> sending_shares(cipher_val.size());
    // std::vector<bool> other_shares(cipher_val.size());
    
    for(size_t i=0; i<cipher_val.size(); i++){
        sending_shares[i] = (char) cipher_val[i].bshares[0];
    }
    runtime.mComm.mPrev.asyncSendCopy(sending_shares.data(), sending_shares.size());
    runtime.mComm.mNext.recv(sending_shares.data(), sending_shares.size());

    std::vector<bool> other_shares(cipher_val.size());
    for(size_t i=0; i<cipher_val.size(); i++){
        other_shares[i] = ((bool) sending_shares[i]) ^ cipher_val[i].bshares[1] ^ cipher_val[i].bshares[0];
    }
    return other_shares;
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