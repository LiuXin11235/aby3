#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>

#ifndef _ABY3_BASICS_H_
#define _ABY3_BASICS_H_

void bool_cipher_lt(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_eq(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_or(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_add(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_and(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void get_permutation(size_t len, std::vector<size_t> &permutation, oc::block &seed);

void get_inverse_permutation(std::vector<size_t> &permutation, std::vector<size_t> &inverse_permutation);

void combine_permutation(std::vector<std::vector<size_t>> &permutation_list, std::vector<size_t> &final_permutation);

template <typename T>
void plain_permutate(std::vector<size_t> &permutation, std::vector<T> &data){
    size_t len = data.size();
    std::vector<T> tmp(len);
    for(size_t i = 0; i < len; i++){
        tmp[permutation[i]] = data[i];
    }
    data = tmp;
}

void get_random_mask(int pIdx, aby3::i64Matrix &res, oc::block &seed);

#endif