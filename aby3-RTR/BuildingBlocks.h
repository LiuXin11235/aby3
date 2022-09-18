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

// template <Decimal D>
// int cipher_mul(u64 pIdx, sf64<D> &a, sf64<D> &b, sf64<D> &res,
//                Sh3Evaluator &eval, Sh3Runtime &runtime);

// template <aby3::Decimal D>
// int cipher_mul(aby3::u64 pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sf64Matrix<D> &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
//     eval.asyncMul(runtime, sharedA, sharedB, res).get();
//     return 0;
// }


// int cipher_mul(aby3::u64 pIdx, aby3::si64Matrix &sharedA, aby3::si64Matrix &sharedB, aby3::si64Matrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
//     eval.asyncMul(runtime, sharedA, sharedB, res).get();
//     return 0;
// }

// synchronized version of fetch_msb.
int fetch_msb(int pIdx, aby3::si64Matrix &diffAB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime, aby3::Sh3Task &task);

// // asynchronized version of fetch_msb.
// int fetch_msb(int pIdx, aby3::si64Matrix &diffAB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime);

// sfixed greater-than.
// template <aby3::Decimal D>
template <aby3::Decimal D>
int cipher_gt(int pIdx, aby3::sf64Matrix<D> &sharedA, aby3::sf64Matrix<D> &sharedB, aby3::sbMatrix &res, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime){
    aby3::sf64Matrix<D> diffF = sharedB - sharedA;
    aby3::si64Matrix diffAB = diffF.i64Cast();
    aby3::Sh3Task task = runtime.noDependencies();
    return fetch_msb(pIdx, diffAB, res, eval, runtime, task);
}

// // sint greater-than.
// int cipher_gt(int pIdx, si64Matrix &sharedA, si64Matrix &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime);

// // sint and plaintext greater-than.
// int cipher_gt(int pIdx, si64Matrix &sharedA, vector<int> &plainB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime);

// // sfixed equal (with wrong outputs).
// template <Decimal D>
// int cipher_eq(u64 pIdx, sf64Matrix<D> &sharedA, sf64Matrix<D> &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime);

// int cipher_eq(u64 pIdx, si64Matrix &intA, si64Matrix &intB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime);