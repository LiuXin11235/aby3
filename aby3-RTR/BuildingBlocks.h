#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>

// setup function.
void basic_setup(aby3::u64 partyIdx, oc::IOService &ios, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval,
           aby3::Sh3Runtime &runtime);

int cipher_mul_seq(int pIdx, const aby3::si64Matrix &sharedA, const aby3::sbMatrix &sharedB, aby3::si64Matrix &res,aby3::Sh3Evaluator &eval, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime &runtime);

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

// // sfixed equal (with wrong outputs).
// template <Decimal D>
// int cipher_eq(u64 pIdx, sf64Matrix<D> &sharedA, sf64Matrix<D> &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime);

int cipher_eq(int pIdx, aby3::si64Matrix &intA, aby3::si64Matrix &intB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);