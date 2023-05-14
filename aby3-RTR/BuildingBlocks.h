#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>

#ifndef _A_H_
#define _A_H_

// setup function.
void distribute_setup(aby3::u64 partyIdx, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
           aby3::Sh3Runtime &runtime);

// local setup funtion, each party is assigned to a thread.
void basic_setup(aby3::u64 partyIdx, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
           aby3::Sh3Runtime &runtime);

int pi_cb_mul(int pIdx, const aby3::i64Matrix &plainA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime &runtime);


int cipher_mul_seq(int pIdx, const aby3::si64Matrix &sharedA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

template <aby3::Decimal D>
int cipher_mul_seq(int pIdx, const aby3::sf64Matrix<D> &sharedA, const aby3::sbMatrix &sharedB, aby3::sf64Matrix<D> &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime){
    return cipher_mul_seq(pIdx, sharedA.i64Cast(), sharedB, res.i64Cast(), eval, enc, runtime);
}

inline int cipher_mul_seq(int pIdx, const aby3::i64Matrix &plainA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime &runtime){
    return pi_cb_mul(pIdx, plainA, sharedB, res, eval, enc, runtime);
}

// synchronized version of fetch_msb.
int fetch_msb(int pIdx, aby3::si64Matrix &diffAB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime, aby3::Sh3Task &task);

// asynchronized version of fetch_msb.
int fetch_msb(int pIdx, aby3::si64Matrix &diffAB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sint greater-than.
int cipher_gt(int pIdx, aby3::si64Matrix &sharedA, aby3::si64Matrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sint and plaintext greater-than.
int cipher_gt(int pIdx, aby3::si64Matrix &sharedA, std::vector<int> &plainB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sfixed greater-than.
template <aby3::Decimal D>
int cipher_gt(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    return cipher_gt(pIdx, sharedA.i64Cast(), sharedB.i64Cast(), res, eval, runtime);
}

// eq implemeneted through two fetch_msb and an Nor gate.
int cipher_eq(int pIdx, aby3::si64Matrix &intA, aby3::si64Matrix &intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);


template <aby3::Decimal D>
int cipher_eq(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    return cipher_eq(pIdx, sharedA.i64Cast(), sharedB.i64Cast(), res, eval, runtime);
}

// eq implemented through int_eq circuit.
int circuit_cipher_eq(int pIdx, aby3::si64Matrix &intA, aby3::si64Matrix &intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

int vector_cipher_eq(int pIdx, std::vector<aby3::si64>& intA, std::vector<int>& intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// fetch_eq_res.
int fetch_eq_res(int pIdx, aby3::sbMatrix& circuitA, aby3::sbMatrix& circuitB, aby3::sbMatrix& res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// initiate functions.
int init_ones(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::si64Matrix &res, int n);

int init_zeros(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime, aby3::si64Matrix &res, int n);

#endif