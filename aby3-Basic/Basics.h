#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>

#ifndef _ABY3_BASICS_H_
#define _ABY3_BASICS_H_

struct boolShare{
    std::array<bool, 2> bshares;

    boolShare(bool share0, bool share1){
        bshares[0] = share0;
        bshares[1] = share1;
    }

    boolShare(){
        bshares[0] = false;
        bshares[1] = false;
    }

    boolShare(bool plain_val, int pIdx){
        switch(pIdx){
            case 0:
                bshares[0] = false;
                bshares[1] = false;
                break;
            case 1:
                bshares[0] = plain_val;
                bshares[1] = false;
                break;
            case 2:
                bshares[0] = false;
                bshares[1] = plain_val;
                break;
            default:
                throw std::runtime_error("boolShare: invalid pIdx");
        }
    }

    aby3::sbMatrix to_matrix(){
        aby3::sbMatrix res(1, 1);
        res.mShares[0](0, 0) = bshares[0];
        res.mShares[1](0, 0) = bshares[1];
        return res;
    }

    void from_matrix(aby3::i64 s1, aby3::i64 s2){
        bshares[0] = s1 & 1; bshares[1] = s2 & 1;
    }
};

struct vecBoolShares{
    std::vector<boolShare> bshares;

    vecBoolShares(std::vector<boolShare>& share_data){
        bshares = share_data;
    }

    vecBoolShares(size_t len){
        bshares.resize(len);
    }

    aby3::sbMatrix to_matrix(){
        size_t len = bshares.size();
        aby3::sbMatrix res(len, 1);
        for(size_t i = 0; i < len; i++){
            res.mShares[0](i, 0) = bshares[i].bshares[0];
            res.mShares[1](i, 0) = bshares[i].bshares[1];
        }
        return res;
    }

    void from_matrix(aby3::sbMatrix &mat){
        size_t len = mat.rows();
        bshares.resize(len);
        for(size_t i = 0; i < len; i++){
            bshares[i].bshares[0] = mat.mShares[0](i, 0) & 1;
            bshares[i].bshares[1] = mat.mShares[1](i, 0) & 1;
        }
    }
};

struct boolIndex{
    std::array<aby3::i64, 2> indexShares;   

    boolIndex(aby3::i64 share0, aby3::i64 share1){
        indexShares[0] = share0;
        indexShares[1] = share1;
    }
    
    boolIndex(){
        indexShares[0] = 0;
        indexShares[1] = 0;
    }

    aby3::sbMatrix to_matrix(){
        aby3::sbMatrix res(1, 1);
        res.mShares[0](0, 0) = indexShares[0];
        res.mShares[1](0, 0) = indexShares[1];
        return res;
    }

    void from_matrix(aby3::sbMatrix &m){
        assert (m.mShares[0].rows() == 1 && m.mShares[0].cols() == 1);
        indexShares[0] = m.mShares[0](0, 0);
        indexShares[1] = m.mShares[1](0, 0);
    }
};

void bool_cipher_lt(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_eq(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_or(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_or(int pIdx, boolShare &sharedA, boolShare &sharedB, boolShare &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_add(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_and(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_and(int pIdx, boolShare &sharedA, boolShare &sharedB, boolShare &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_cipher_not(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &res);

void bool_cipher_not(int pIdx, std::vector<boolShare> &sharedA, std::vector<boolShare> &res);

void bool_cipher_not(int pIdx, boolShare &sharedA, boolShare &res);

void bool_cipher_dot(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_get_first_zero_mask(int pIdx, std::vector<boolShare>& inputA, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

void bool_init_false(int pIdx, aby3::sbMatrix &res);

void bool_init_true(int pIdx, aby3::sbMatrix &res);

void bool_init_false(int pIdx, boolShare &res);

void bool_init_true(int pIdx, boolShare &res);

void bool_shift_and_left(int pIdx, aby3::sbMatrix &sharedA, size_t shift_len, aby3::sbMatrix &res_shift, aby3::sbMatrix &res_left);

aby3::i64Matrix back2plain(int pIdx, aby3::sbMatrix &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

bool back2plain(int pIdx, boolShare &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

aby3::i64 back2plain(int pIdx, boolIndex &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

std::vector<bool> back2plain(int pIdx, std::vector<boolShare> &cipher_val, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

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