#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <map>

#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>

#include "aby3/sh3/Sh3Encryptor.h"
#include "aby3/sh3/Sh3Evaluator.h"
#include "aby3/sh3/Sh3FixedPoint.h"
#include "aby3/sh3/Sh3Runtime.h"
#include "aby3/sh3/Sh3Types.h"
#include "aby3/sh3/Sh3BinaryEvaluator.h"
#include "aby3/Circuit/CircuitLibrary.h"

using namespace std;
using namespace oc;
using namespace aby3;

template <Decimal D>
int cipher_mul(u64 pIdx, sf64<D> &a, sf64<D> &b, sf64<D> &res,
               Sh3Evaluator &eval, Sh3Runtime &runtime) {
  eval.asyncMul(runtime, a, b, res).get();
  return 0;
}

template <Decimal D>
int cipher_mul(u64 pIdx, sf64Matrix<D> &sharedA, sf64Matrix<D> &sharedB, sf64Matrix<D> &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  eval.asyncMul(runtime, sharedA, sharedB, res).get();
  return 0;
}

int fetch_msb(int pIdx, si64Matrix &diffAB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  // 1. set the binary circuits.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(diffAB.size(), 64);
  circuitInput1.resize(diffAB.size(), 64);

  // 2. let party0 adds up two shares x0 and x1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<diffAB.size(); j++){
        circuitInput0.mShares[0](j) = diffAB.mShares[0](j)+ diffAB.mShares[1](j);
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;
    }
    case 1: {
      circuitInput1.mShares[0].resize(diffAB.rows(), diffAB.cols());
      circuitInput1.mShares[1].setZero();
      memcpy(circuitInput1.mShares[0].data(), diffAB.mShares[0].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].resize(diffAB.rows(), diffAB.cols());
      memcpy(circuitInput1.mShares[1].data(), diffAB.mShares[1].data(), diffAB.size() * sizeof(i64));

      circuitInput0.mShares[0].setZero();
    }
  }

  // 3. binary secrets reshare.
  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(), circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(), circuitInput0.mShares[1].size());
  fu.get();

  // 4. extract the msb of the difference.
  auto cir = lib.int_comp_helper(sizeof(i64)*8);
  binEng.setCir(cir, diffAB.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);

  binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
    res.resize(diffAB.size(), 1);
    binEng.getOutput(0, res);
  }).get();

  return 0;
}


template <Decimal D>
int cipher_gt(int pIdx, sf64Matrix<D> &sharedA, sf64Matrix<D> &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // comparision between sharedA and sharedB, diffAB = sharedA - sharedB, sub by sint type.
  sf64Matrix<D> diffF = sharedB - sharedA;
  si64Matrix diffAB = diffF.i64Cast();

  return fetch_msb(pIdx, diffAB, res, eval, runtime);
}

int cipher_gt(int pIdx, si64Matrix &sharedA, si64Matrix &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // comparision between sharedA and sharedB, diffAB = sharedA - sharedB, sub by sint type.
  si64Matrix diffAB = sharedB - sharedA;
  return fetch_msb(pIdx, diffAB, res, eval, runtime);
}

int cipher_gt(int pIdx, si64Matrix &sharedA, vector<int> &plainB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  
  // 1. set the binary circuits.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(sharedA.size(), 64);
  circuitInput1.resize(sharedA.size(), 64);

  // 2. let party0 adds up two shares x0 and x1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<sharedA.size(); j++){
        circuitInput0.mShares[0](j) = plainB[j] - (sharedA.mShares[0](j)+ sharedA.mShares[1](j));
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;
    }
    case 1: {
      circuitInput1.mShares[0].resize(sharedA.rows(), sharedA.cols());
      circuitInput1.mShares[1].setZero();
      memcpy(circuitInput1.mShares[0].data(), sharedA.mShares[0].data(), sharedA.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].resize(sharedA.rows(), sharedA.cols());
      memcpy(circuitInput1.mShares[1].data(), sharedA.mShares[1].data(), sharedA.size() * sizeof(i64));

      circuitInput0.mShares[0].setZero();
    }
  }

  // 3. binary secrets reshare.
  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(), circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(), circuitInput0.mShares[1].size());
  fu.get();

  // 4. extract the msb of the difference.
  auto cir = lib.int_comp_helper(sizeof(i64)*8);
  binEng.setCir(cir, sharedA.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);

  binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
    res.resize(sharedA.size(), 1);
    binEng.getOutput(0, res);
  }).get();

  return 0;
}


template <Decimal D>
int cipher_eq(u64 pIdx, sf64Matrix<D> &sharedA, sf64Matrix<D> &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  si64Matrix intA = sharedA.i64Cast();
  si64Matrix intB = sharedB.i64Cast();
  
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;

  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(intA.size(), 64);
  circuitInput1.resize(intB.size(), 64);
  memcpy(circuitInput0.mShares[0].data(), intA.mShares[0].data(), intA.size() * sizeof(i64));
  memcpy(circuitInput0.mShares[1].data(), intA.mShares[1].data(), intA.size() * sizeof(i64));
  memcpy(circuitInput1.mShares[0].data(), intB.mShares[0].data(), intB.size() * sizeof(i64));
  memcpy(circuitInput1.mShares[1].data(), intB.mShares[1].data(), intB.size() * sizeof(i64));

  // 3. generate the eq circuit and evaluate.
  auto cir = lib.int_eq(sizeof(i64)*8);
  binEng.setCir(cir, intA.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);
  binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
    res.resize(intA.size(), 1);
    binEng.getOutput(0, res);
  }).get();

  return 0;
}


int cipher_eq(u64 pIdx, si64Matrix &intA, si64Matrix &intB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // 1. set the correspondong boolean circuit.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;

  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(intA.size(), 64);
  circuitInput1.resize(intB.size(), 64);
  memcpy(circuitInput0.mShares[0].data(), intA.mShares[0].data(), intA.size() * sizeof(i64));
  memcpy(circuitInput0.mShares[1].data(), intA.mShares[1].data(), intA.size() * sizeof(i64));
  memcpy(circuitInput1.mShares[0].data(), intB.mShares[0].data(), intB.size() * sizeof(i64));
  memcpy(circuitInput1.mShares[1].data(), intB.mShares[1].data(), intB.size() * sizeof(i64));

  // 3. generate the eq circuit and evaluate.
  auto cir = lib.int_eq(sizeof(i64)*8);
  binEng.setCir(cir, intA.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);
  binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
    res.resize(intA.size(), 1);
    binEng.getOutput(0, res);
  }).get();

  return 0;
}


template <Decimal D>
int cipher_index(u64 pIdx, sf64Matrix<D> &sharedM, si64Matrix &cipherIndex, sf64Matrix<D> &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  return 0;
}

// cipher_argsort_offset, try to implement the blockwise computation.
template <Decimal D>
int cipher_argsort_offset(int pIdx, sf64Matrix<D> &sharedM, si64Matrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime, Sh3Encryptor &enc, int offsetLeft, int offsetRight){
  // 1. compute the diff between target elements from offsetLeft to offsetRight.
  u64 block_length = offsetRight - offsetLeft;
  u64 n = sharedM.rows();
  si64Matrix diff(block_length, 1);
  si64Matrix isharedM = sharedM.i64Cast();

  u64 start_i = offsetLeft / n, start_j = offsetLeft % n;
  u64 end_i = offsetRight / n, end_j = offsetRight % n;
  u64 p = 0, i = start_i, j = start_j;

  while (i < end_i) {
    for (j; j < n; j++) {
      diff.mShares[0](p, 0) = isharedM.mShares[0](j, 0) - isharedM.mShares[0](i, 0);
      diff.mShares[1](p, 0) = isharedM.mShares[1](j, 0) - isharedM.mShares[1](i, 0);
      p++;
    }
    j = 0;
    i++;
  }
  for (j = 0; j < end_j; j++) {
    diff.mShares[0](p, 0) = isharedM.mShares[0](j, 0) - isharedM.mShares[0](end_i, 0);
    diff.mShares[1](p, 0) = isharedM.mShares[1](j, 0) - isharedM.mShares[1](end_i, 0);
    p++;
  }
  sbMatrix pairwise_comp;
  fetch_msb(pIdx, diff, pairwise_comp, eval, runtime);
  cout << "succeed fetch msb" << endl;
  // retrive the pairwise comparision result and convert to secret ints.
  i64Matrix ones(block_length, 1);
  for(int i=0; i<block_length; i++) ones(i, 0) = 1;
  si64Matrix ipairwise(block_length, 1);
  if(pIdx == 0){
    enc.localIntMatrix(runtime, ones, ipairwise).get();
  }else{
    enc.remoteIntMatrix(runtime, ipairwise).get();
  }
  eval.asyncMul(runtime, ipairwise, pairwise_comp, ipairwise).get();
  // reduce the result to the global res.
  i = start_i, j = start_j;
  p = 0;
  while (i < end_i) {
    for (j; j < n; j++) {
      res.mShares[0](i, 0) += ipairwise.mShares[0](p, 0);
      res.mShares[1](i, 0) += ipairwise.mShares[1](p, 0);
      p++;
    }
    j = 0;
    i++;
  }
  for (j = 0; j < end_j; j++) {
    res.mShares[0](i, 0) += ipairwise.mShares[0](p, 0);
    res.mShares[1](i, 0) += ipairwise.mShares[1](p, 0);
    p++;
  }

  return 0;
}

template <Decimal D>
int cipher_argsort(int pIdx, sf64Matrix<D> &sharedM, si64Matrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime, Sh3Encryptor &enc){

  // 1. expand sharedM to sharedM_expansion, repeat the n elements n times.
  u64 n = sharedM.rows();
  u64 expand_length = n * n;
  sf64Matrix<D> sharedM_expand(expand_length, 1), sharedM_expandT(expand_length, 1);

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      sharedM_expand(i * n + j, 0, sharedM(i, 0));
      sharedM_expandT(i * n + j, 0, sharedM(j, 0));
    }
  }

  // 2. compare sharedM_expand and sharedM_expandT to get the pairwise comparision result.
  sbMatrix comp_res;
  if(cipher_gt(pIdx, sharedM_expand, sharedM_expandT, comp_res, eval, runtime) != 0){
    throw runtime_error("something happened in cipher_gt.  "LOCATION);
  }

  // 3. transfer the sbMatrix comp_res to siMatrix and compute the final result.
  i64Matrix ones(expand_length, 1);
  i64Matrix zeros(n, 1);
  for(int i=0; i<expand_length; i++) ones(i, 0) = 1;
  for(int i=0; i<n; i++) zeros(i, 0) = 0;
  si64Matrix icomp_res(expand_length, 1);
  res.resize(n, 1);

  if(pIdx == 0){
    enc.localIntMatrix(runtime, ones, icomp_res).get();
    enc.localIntMatrix(runtime, zeros, res).get();
  }
  else{
    enc.remoteIntMatrix(runtime, icomp_res).get();
    enc.remoteIntMatrix(runtime, res).get();
  }
  eval.asyncMul(runtime, icomp_res, comp_res, icomp_res).get();

  // 4. reduce to the final result;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      res.mShares[0](i, 0) += icomp_res.mShares[0](i * n + j, 0);
      res.mShares[1](i, 0) += icomp_res.mShares[1](i * n + j, 0);
    }
  }

  return 0;
}


template <Decimal D>
int rtr_cipher_argsort(int pIdx, sf64Matrix<D> &sharedM, si64Matrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime, Sh3Encryptor &enc){

  u64 n = sharedM.rows();
  u64 expand_length = n * n;
  int offsetLeft = 0, offsetRight = expand_length;

  // initialize res to all zeros matrix.
  i64Matrix zeros(n, 1);
  for(int i=0; i<n; i++) zeros(i, 0) = 0;
  res.resize(n, 1);
  if(pIdx == 0){
    enc.localIntMatrix(runtime, zeros, res).get();
  }
  else{
    enc.remoteIntMatrix(runtime, res).get();
  }
  
  return cipher_argsort_offset(pIdx, sharedM, res, eval, runtime, enc, offsetLeft, offsetRight);
}


// // cipher_argsort_offset, try to implement the blockwise computation.
// template <Decimal D>
// int cipher_argsort_offset(int pIdx, sf64Matrix<D> &sharedM, si64Matrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime, Sh3Encryptor &enc, int offsetLeft, int offsetRight){
//   // 1. compute the diff between target elements from offsetLeft to offsetRight.
//   u64 block_length = offsetRight - offsetLeft;
//   u64 n = sharedM.rows();
//   si64Matrix diff(block_length, 1);
//   si64Matrix isharedM = sharedM.i64Cast();

//   u64 start_i = offsetLeft / n, start_j = offsetLeft % n;
//   u64 end_i = offsetRight / n, end_j = offsetRight % n;
//   u64 p = 0, i = start_i, j = start_j;

//   while (i < end_i) {
//     for (j; j < n; j++) {
//       diff.mShares[0](p, 0) = isharedM.mShares[0](i, 0) - isharedM.mShares[0](j, 0);
//       diff.mShares[1](p, 0) = isharedM.mShares[1](i, 0) - isharedM.mShares[1](j, 0);
//       p++;
//     }
//     j = 0;
//     i++;
//   }
//   for (j = 0; j < end_j; j++) {
//     diff.mShares[0](p, 0) = isharedM.mShares[0](end_i, 0) - isharedM.mShares[0](j, 0);
//     diff.mShares[1](p, 0) = isharedM.mShares[1](end_i, 0) - isharedM.mShares[1](j, 0);
//     p++;
//   }
//   sbMatrix pairwise_comp;
//   fetch_msb(pIdx, diff, pairwise_comp, eval, runtime);

//   // retrive the pairwise comparision result and convert to secret ints.
//   i64Matrix ones(block_length, 1);
//   for(int i=0; i<block_length; i++) ones(i, 0) = 1;
//   si64Matrix ipairwise(block_length, 1);
//   if(pIdx == 0){
//     enc.localIntMatrix(runtime, ones, ipairwise).get();
//   }else{
//     enc.remoteIntMatrix(runtime, ipairwise).get();
//   }
  
//   // reduce the result to the global res.
//   i = start_i, j = start_j;
//   p = 0;
//   while (i < end_i) {
//     for (j; j < n; j++) {
//       res.mShares[0](i, 0) += ipairwise.mShares[0](i * n + j, 0);
//       res.mShares[1](i, 0) += ipairwise.mShares[1](i * n + j, 0);
//       // res_x[i] += pairwise_comp.x[p];
//       // res_x_[i] += pairwise_comp.x_[p];
//       p++;
//     }
//     j = 0;
//     i++;
//   }
//   for (j = 0; j < end_j; j++) {
//     res.mShares[0](i, 0) += ipairwise.mShares[0](i * n + j, 0);
//     res.mShares[1](i, 0) += ipairwise.mShares[1](i * n + j, 0);
//     // res_x[end_i] += pairwise_comp.x[p];
//     // res_x_[end_i] += pairwise_comp.x_[p];
//     p++;
//   }

//   return 0;
// }


void basic_setup(u64 partyIdx, IOService &ios, Sh3Encryptor &enc, Sh3Evaluator &eval,
           Sh3Runtime &runtime) {
  CommPkg comm;
  switch (partyIdx) {
    case 0:
      comm.mNext = Session(ios, "127.0.0.1:1313", SessionMode::Server, "01")
                       .addChannel();
      comm.mPrev = Session(ios, "127.0.0.1:1314", SessionMode::Server, "02")
                       .addChannel();
      break;
    case 1:
      comm.mNext = Session(ios, "127.0.0.1:1315", SessionMode::Server, "12")
                       .addChannel();
      comm.mPrev = Session(ios, "127.0.0.1:1313", SessionMode::Client, "01")
                       .addChannel();
      break;
    default:
      comm.mNext = Session(ios, "127.0.0.1:1314", SessionMode::Client, "02")
                       .addChannel();
      comm.mPrev = Session(ios, "127.0.0.1:1315", SessionMode::Client, "12")
                       .addChannel();
      break;
  }
    // Establishes some shared randomness needed for the later protocols
    enc.init(partyIdx, comm, sysRandomSeed());

    // Establishes some shared randomness needed for the later protocols
    eval.init(partyIdx, comm, sysRandomSeed());

    // Copies the Channels and will use them for later protcols.
    runtime.init(partyIdx, comm);
}


int test_mul(CLP &cmd) {
  IOService ios;
  f64<D8> a = 2.2, b = 3.3;

  vector<thread> thrds;
  for (u64 i = 0; i < 3; i++) {
    thrds.emplace_back(thread([i, a, b, &ios]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup(i, ios, enc, eval, runtime); // same setup function in the tutorial.
      sf64<D8> ea, eb, res;
      f64<D8> pres;
      if (i == 0) {
        enc.localFixed(runtime, a, ea).get();
        enc.localFixed(runtime, b, eb).get();
      } else {
        enc.remoteFixed(runtime, ea).get();
        enc.remoteFixed(runtime, eb).get();
      }

      cipher_mul(i, ea, eb, res, eval, runtime);
      enc.revealAll(runtime, res, pres).get();
      if(i == 0)
        cout << "final res = " << (double) pres << endl;
    }));
  }
  for (auto &t : thrds)
    t.join();
  return 0;
}

int test_gt(CLP& cmd){
  IOService ios;
  u64 rows = 10, cols = 1;
  f64Matrix<D8> plainA(rows, cols);
  f64Matrix<D8> plainB(rows, cols);

  for(u64 j=0; j<rows; j++){
    for (u64 i = 0; i<cols; i++){
      plainA(j, i) = (double) j;
      plainB(j, i) = (double) j + 1;
    }
  }

  vector<thread> thrds;
  for(int i=0; i<3; i++){
    thrds.emplace_back(thread([i, plainA, plainB, &ios]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup((u64)i, ios, enc, eval, runtime); // same setup function in the tutorial.

      sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
      sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());

      if(i == 0){
        enc.localFixedMatrix(runtime, plainA, sharedA).get();
        enc.localFixedMatrix(runtime, plainB, sharedB).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
      }

      sbMatrix lt_res;
      cipher_gt(i, sharedB, sharedA, lt_res, eval, runtime);
      
      // si64Matrix larger_res;
      sf64Matrix<D8> res;
      res.resize(sharedA.rows(), sharedA.cols());
      // larger_res.resize(sharedA.rows(), sharedA.cols());
      eval.asyncMul(runtime, sharedA.i64Cast(), lt_res, res.i64Cast()).get();

      
      f64Matrix<D8> pres;
      enc.revealAll(runtime, res, pres).get();
      
      if(i == 0){
        cout << "expected res = [0.0, 1.0, 2.0, ..., 9.0]" << endl;
        cout << "comp and extract: ";
        for(int j=0; j<sharedA.rows(); j++) cout << pres(j, 0) << ", ";
        cout << endl;
        // cout << "comp and extract = " << pres << endl;
      }
    }));
  }

  for (auto &t : thrds)
    t.join();
  return 0;
}

int test_eq(CLP& cmd){
  IOService ios;
  u64 rows = 10, cols = 1;

  // test for f64 datatype.
  f64Matrix<D8> plainA(rows, cols);
  f64Matrix<D8> plainB(rows, cols);
  f64Matrix<D8> ones(rows, cols);

  for(u64 j=0; j<rows; j++){
    for (u64 i = 0; i<cols; i++){
      plainA(j, i) = (double) j;
      ones(j, i) = 1;
      if (j%2 == 0) plainB(j, i) = (double) j + 1;
      else plainB(j, i) = (double) j;
    }
  } 

  // test for i64 datatype.
  i64Matrix iplainA(rows, cols);
  i64Matrix iplainB(rows, cols);
  i64Matrix iones(rows, cols);

  for(u64 j=0; j<rows; j++){
    for (u64 i = 0; i<cols; i++){
      iplainA(j, i) = j;
      iones(j, i) = 1;
      // iplainB(j, i) = j;
      if (j%2 == 0) iplainB(j, i) = j + 1;
      else iplainB(j, i) = j;
    }
  }

  vector<thread> thrds;
  for(u64 i=0; i<3; i++){
    thrds.emplace_back(thread([i, plainA, plainB, ones, iplainA, iplainB, iones, &ios]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup(i, ios, enc, eval, runtime); // same setup function in the tutorial.

      sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
      sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());
      sf64Matrix<D8> sharedOnes(plainA.rows(), plainA.cols());

      si64Matrix isharedA(plainA.rows(), plainA.cols());
      si64Matrix isharedB(plainA.rows(), plainA.cols());
      si64Matrix isharedOnes(plainA.rows(), plainA.cols());

      if(i == 0){
        enc.localFixedMatrix(runtime, plainA, sharedA).get();
        enc.localFixedMatrix(runtime, plainB, sharedB).get();
        enc.localFixedMatrix(runtime, ones, sharedOnes).get();

        enc.localIntMatrix(runtime, iplainA, isharedA).get();
        enc.localIntMatrix(runtime, iplainB, isharedB).get();
        enc.localIntMatrix(runtime, iplainA, isharedOnes).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
        enc.remoteFixedMatrix(runtime, sharedOnes).get();

        enc.remoteIntMatrix(runtime, isharedA).get();
        enc.remoteIntMatrix(runtime, isharedB).get();
        enc.remoteIntMatrix(runtime, isharedOnes).get();
      }

      sbMatrix eq_res;
      cipher_eq(i, sharedA, sharedB, eq_res, eval, runtime);
      sf64Matrix<D8> res;
      res.resize(sharedA.rows(), sharedA.cols()); 
      eval.asyncMul(runtime, sharedOnes.i64Cast(), eq_res, res.i64Cast()).get();
      f64Matrix<D8> pres;
      enc.revealAll(runtime, res, pres).get();
      
      if(i == 0){
        cout << "expected res = [1.0, 0.0, 1.0, ...]" << endl;
        cout << "equal test = " << pres << endl;
      }

      sbMatrix ieq_res;
      cipher_eq(i, isharedA, isharedB, ieq_res, eval, runtime);

      si64Matrix ires;
      ires.resize(isharedA.rows(), isharedA.cols());
      eval.asyncMul(runtime, isharedOnes, ieq_res, ires).get();
      i64Matrix ipres;
      enc.revealAll(runtime, ires, ipres).get();
      if(i == 0){
        cout << "test si equal" << endl;
        cout << "expected res = [1, 0, 1, ...]" << endl;
        cout << "equal test = " << ipres << endl;
      }

    }));
  }

  for (auto &t : thrds)
    t.join();
  return 0;
}

int test_argsort(CLP& cmd){
  IOService ios;

  // generate the test data.
  u64 rows = 5, cols = 1;
  double test[rows] = {3, 5, 7, 6, 9};
  f64Matrix<D8> plainTest(rows, cols);
  for(int i=0; i<rows; i++){
    plainTest(i, 0) = test[i];
  }

  // start the three parties computation.
  vector<thread> thrds;
  for(u64 i=0; i<3; i++){
    thrds.emplace_back(thread([i, &ios, plainTest, rows, cols]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup(i, ios, enc, eval, runtime);

      // generate the cipher test data.
      sf64Matrix<D8> sharedM(rows, cols);
      if(i == 0){
        enc.localFixedMatrix(runtime, plainTest, sharedM).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedM).get();
      }
      
      si64Matrix argres;
      // cipher_argsort(i, sharedM, argres, eval, runtime, enc);
      rtr_cipher_argsort(i, sharedM, argres, eval, runtime, enc);

      // reveal the result and check, expected res: [0, 1, 3, 2, 4].
      i64Matrix plain_argres;
      enc.revealAll(runtime, argres, plain_argres).get();
      if(i == 0){
        cout << "Expected result: [0 1 3 2 4]" << endl;
        // cout << plain_argres << endl;
        for(int j=0; j<rows; j++){
          cout << plain_argres(j, 0) << " ";
        }
        cout << endl;
      }
    }));
  }
  for (auto &t : thrds)
    t.join();
  return 0;
}

int basic_performance(CLP& cmd, int n, int repeats, map<string, vector<double>>& dict){
  IOService ios;
  u64 rows = n, cols = 1;
  f64Matrix<D8> plainA(rows, cols);
  f64Matrix<D8> plainB(rows, cols);

  for(u64 j=0; j<rows; j++){
    for (u64 i = 0; i<cols; i++){
      plainA(j, i) = (double) j;
      plainB(j, i) = (double) j + 1;
    }
  }

  vector<double> time_mul_array(3);
  vector<double> time_gt_array(3);

  vector<thread> thrds;
  for(u64 i=0; i<3; i++){
    thrds.emplace_back(thread([i, plainA, plainB, repeats, &time_mul_array, &time_gt_array, &ios]() {
      
      // setup the environment.
      clock_t start, end;
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup(i, ios, enc, eval, runtime);

      // generate the test matrix.
      sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
      sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());

      if(i == 0){
        enc.localFixedMatrix(runtime, plainA, sharedA).get();
        enc.localFixedMatrix(runtime, plainB, sharedB).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
      }

      // test the performance of basic ops.
      sbMatrix lt_res;
      start = clock();
      for(int k=0; k<repeats; k++)
        cipher_gt(i, sharedB, sharedA, lt_res, eval, runtime);
      end = clock();

      double time_gt = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
      time_gt_array[i] = time_gt;
      // cout << "time_gt_array: " << time_gt_array[i] << endl;

      sf64Matrix<D8> mul_res(plainA.rows(), plainA.cols());
      start = clock();
      for(int k=0; k<repeats; k++)
        cipher_mul(i, sharedA, sharedB, mul_res, eval, runtime);
      end = clock();

      double time_mul = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
      time_mul_array[i] = time_mul;

    }));
  }

  for (auto &t : thrds)
    t.join();
  
  // map<string, double(*)[3]> dict;
  dict["mul"] = time_mul_array;
  dict["gt"] = time_gt_array;

  // for(int i=0; i<3; i++) cout << dict["mul"][i] << endl;

  return 0;
}