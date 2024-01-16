#include "Basics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <random>

#include <aby3/sh3/Sh3BinaryEvaluator.h>
#include <aby3/Circuit/CircuitLibrary.h>

#include "../aby3-RTR/debug.h"

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